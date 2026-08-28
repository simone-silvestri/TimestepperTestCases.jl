using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using NumericalEarth
using Dates: DateTime
using Printf

using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces:
    averaging_shape_function,
    LowDissipationAveragingKernel,
    SymmetricTrigAveragingKernel,
    WideTrig74AveragingKernel,
    WideTrig2AveragingKernel,
    OptimizedAsymmetricAveragingKernel,
    ForwardBackwardScheme,
    RungeKutta3Scheme

using Oceananigans.BuoyancyFormulations: LinearEquationOfState
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation

using Downloads: Downloads
using NumericalEarth.DataWrangling: metadata_path

const ARTIFACTS_BASE_URL = "https://github.com/NumericalEarth/NumericalEarthArtifacts/releases/download/data-v1/"

function download_from_artifacts(filepath)
    if isfile(filepath)
        return
    end

    filename = basename(filepath)
    fallback_url = ARTIFACTS_BASE_URL * filename
    @info "Downloading $filename from NumericalEarthArtifacts fallback..."
    mktemp(dirname(filepath)) do tmppath, tmpio
        close(tmpio)
        Downloads.download(fallback_url, tmppath)
        mv(tmppath, filepath; force=true)
    end
end

NumericalEarth.EarthSystemModels.reference_density(::LinearEquationOfState) = 1026
NumericalEarth.EarthSystemModels.heat_capacity(::LinearEquationOfState) = 3992

near_global_kernel(name::Symbol)    = near_global_kernel(Val(name))
near_global_kernel(::Val{:SM05})    = averaging_shape_function
near_global_kernel(::Val{:mu2})     = LowDissipationAveragingKernel()
near_global_kernel(::Val{:trig})    = SymmetricTrigAveragingKernel()
near_global_kernel(::Val{:trig74})  = WideTrig74AveragingKernel()
near_global_kernel(::Val{:trig2})   = WideTrig2AveragingKernel()
near_global_kernel(::Val{:optasym}) = OptimizedAsymmetricAveragingKernel()

function near_global_grid(arch = CPU();
                          Nx = 1440,
                          Ny = 600,
                          Nz = 50,
                          depth = 5000meters,
                          latitude  = (-75, 75),
                          longitude = (0, 360),
                          minimum_depth = 15meters,
                          interpolation_passes = 5,
                          major_basins = 1)

    z = ExponentialDiscretization(Nz, -depth, 0, mutable=true)

    grid = LatitudeLongitudeGrid(arch;
                                 size = (Nx, Ny, Nz),
                                 halo = (7, 7, 7),
                                 z, latitude, longitude)

    bottom_height = regrid_bathymetry(grid; minimum_depth, interpolation_passes, major_basins)

    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)
end

# The near-global time steps are empirical rather than derived. The theoretical criterion used by the idealized
# cases, Δt = 0.7 θ★/(c₁ k), is not well posed here: both the first baroclinic speed and the grid spacing vary
# over the domain (Δx shrinks as cos φ, from ≈ 27.8 km at the equator to ≈ 7.2 km at 75°), so the binding
# combination c₁ k is a maximum over the globe and needs the actual stratification field to evaluate. The
# three-stage value below is measured; the other schemes are scaled from it by the ratio of the imaginary-axis
# limits, which is the part of the criterion that does transfer.
near_global_timestep(::Val{:SplitRungeKutta3}) = 20minutes

# Every other scheme is derived from that single reference by the ratio of the imaginary-axis limits
# θ★(scheme)/θ★(WRK3), exactly as in the idealized cases: AB2 gives up a factor 0.29 and MRK4 earns 1.591.
near_global_timestep(::Val{scheme}) where scheme =
    baroclinic_timestep(scheme, near_global_timestep(Val(:SplitRungeKutta3)))

"""
    near_global_stability_parameters(grid)

Depth and horizontal spacings that set the barotropic Courant number of this case, read off the grid rather
than stated as numbers: the spacing varies over the domain, so there is no single value to quote.

`Δx` is the zonal spacing at the poleward edge, ≈ 7.2 km at 75° against ≈ 27.8 km at the equator, `Δy` the
meridional one, and `H` the deepest column, so `c₀ k` is evaluated where it binds.
"""
function near_global_stability_parameters(grid)
    return (; H = grid.Lz, Δx = minimum_xspacing(grid), Δy = minimum_yspacing(grid))
end

"""
    near_global_substeps(barotropic_scheme, grid, scheme; averaging_kernel, Δt)

Barotropic substep count for this case, fixed so that the substep Courant number `c₀ k Δτ` sits at 70% of the
limit of `barotropic_scheme` -- `√3` for the three-stage Runge-Kutta substep, `1` for forward-backward -- as in
every other case. The baroclinic step is empirical here, but the substep count that goes with it is derived, so
the barotropic Courant number stays a stated property of the run.

Contrarily to the idealized cases, `k` is the staggered wavenumber `2√(Δx⁻² + Δy⁻²)` of
[`staggered_wavenumber`](@ref) and not the spectral `π/Δx`. The grid is anisotropic where it is tightest,
`Δy/Δx ≈ 3.9` at 75°, so the two conventions no longer agree to the 11% they agree to on an isotropic grid:
here the spectral one overstates the Courant number of the two-point barotropic operator by 1.52.

The grid is a positional argument because the configuration is described by the grid itself, `near_global_grid`
carrying resolution and depth as keyword arguments.
"""
function near_global_substeps(barotropic_scheme, grid, scheme = :SplitRungeKutta3;
                              averaging_kernel = OptimizedAsymmetricAveragingKernel(),
                              Δt = near_global_timestep(Val(scheme)))
    p = near_global_stability_parameters(grid)
    return barotropic_substeps(barotropic_scheme; p.H, p.Δx, Δt, averaging_kernel,
                               wavenumber = staggered_wavenumber(p.Δx, p.Δy))
end

"""
    near_global_discretizations()

The cases of Table 1 run on the near-global configuration: [`discretizations`](@ref) without `WRK3-UP`.

The near-global case fixes tracer advection at the vertically implicit WENO7 of [`near_global`](@ref), whereas
`WRK3-UP` carries the explicitly discretized third-order upwind it shares with the idealized cases. Run here it
would depart from `WRK3-SE` in the vertical time discretization as well as in the spatial reconstruction, so it
no longer isolates the spatial contribution and is left to the idealized cases.
"""
near_global_discretizations() = filter(d -> d.label != "WRK3-UP", discretizations())

"""
    near_global_closure(FT = Oceananigans.defaults.FloatType)

The CATKE closure of this case: the vertically implicit `CATKEVerticalDiffusivity` with the bottom distance
coefficient of the shear mixing length set to `Cᵇ = 0.01`, against the 0.28 of `CATKEMixingLength`.
"""
near_global_closure(FT = Oceananigans.defaults.FloatType) =
    CATKEVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), FT;
                             mixing_length = CATKEMixingLength(Cᵇ = 0.01),
                             turbulent_kinetic_energy_equation = CATKEEquation(Cᵂϵ = 1.0))

function near_global(timestepper::Symbol = :SplitRungeKutta3;
                     arch = CPU(),
                     grid = near_global_grid(arch),
                     free_surface = nothing,
                     filter::Symbol = :optasym,
                     averaging_kernel = near_global_kernel(filter),
                     barotropic_timestepper = ForwardBackwardScheme(),
                     slow_forcing = FrozenSlowForcing(),
                     tracer_advection = nothing,
                     Δt = near_global_timestep(Val(timestepper)),
                     cold_start_Δt = Δt / 3,
                     cold_start_duration = 60days,
                     stop_time = 720days,
                     dissipation = true,
                     closure = near_global_closure(),
                     bottom_drag_coefficient = 0.003,
                     # uᵦ = 0.1 m s⁻¹ stands for the barotropic tide, which this configuration does not force
                     bottom_drag_background_velocity = 0.1,
                     equation_of_state = LinearEquationOfState(),
                     label = nothing,
                     init_date = DateTime(1993, 1, 1),
                     progress_interval = TimeInterval(5days),
                     surface_output_interval = TimeInterval(1days),
                     dissipation_output_interval = AveragedTimeInterval(30days))

    if free_surface === nothing
        # Fixed substep count, set from the substep integrator's own stability limit rather than from a cfl
        # target, so that the barotropic Courant number is a stated property of the run.
        substeps = near_global_substeps(barotropic_timestepper, grid; averaging_kernel, Δt)
        free_surface = SplitExplicitFreeSurface(grid; substeps, averaging_kernel,
                                                timestepper = barotropic_timestepper,
                                                slow_forcing)
    end
    
    time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5)
    momentum_advection = WENOVectorInvariant(; time_discretization)
    tracer_advection = something(tracer_advection, WENO(order=7; minimum_buffer_upwind_order=3, time_discretization))

    ocean = ocean_simulation(grid; free_surface, timestepper, Δt = cold_start_Δt, stop_time = cold_start_duration,
                             closure, equation_of_state, momentum_advection, tracer_advection,
                             bottom_drag_coefficient, bottom_drag_background_velocity)

    Tmetadata = Metadatum(:temperature, dataset=ECCO2Daily(), date=init_date)
    Smetadata = Metadatum(:salinity,    dataset=ECCO2Daily(), date=init_date)
    
    download_from_artifacts(metadata_path(Tmetadata))
    download_from_artifacts(metadata_path(Smetadata))

    set!(ocean.model, T=Tmetadata, S=Smetadata)

    parent(ocean.model.tracers.T) .= max.(parent(ocean.model.tracers.T), -1.8)

    if dissipation
        VFCC = Oceananigans.AbstractOperations.grid_metric_operation((Face,   Center, Center), Oceananigans.Operators.volume, grid)
        VCFC = Oceananigans.AbstractOperations.grid_metric_operation((Center, Face,   Center), Oceananigans.Operators.volume, grid)
        VCCF = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Face),   Oceananigans.Operators.volume, grid)

        b   = Oceananigans.Models.buoyancy_operation(ocean.model)
        Gbx = ∂x(b)^2 * VFCC
        Gby = ∂y(b)^2 * VCFC
        Gbz = ∂z(b)^2 * VCCF

        ϵb = BuoyancyVarianceDissipation(grid)
        add_callback!(ocean, ϵb, IterationInterval(1))

        diss = merge(BuoyancyVarianceDissipationComputations.flatten_dissipation_fields(ϵb), (; Gbx, Gby, Gbz))
        ocean.output_writers[:dissipation] = JLD2Writer(ocean.model, diss;
                                                        schedule = dissipation_output_interval,
                                                        filename = near_global_filename(label, timestepper, filter, free_surface) * "_dissipation",
                                                        overwrite_existing = true)
    end

    surface = merge(ocean.model.tracers, ocean.model.velocities)
    ocean.output_writers[:surface] = JLD2Writer(ocean.model, surface;
                                                schedule = surface_output_interval,
                                                indices = (:, :, grid.Nz),
                                                filename = near_global_filename(label, timestepper, filter, free_surface) * "_surface",
                                                with_halos = true,
                                                overwrite_existing = true,
                                                array_type = Array{Float32})

    ocean.output_writers[:average] = JLD2Writer(ocean.model, surface;
                                                schedule = dissipation_output_interval,
                                                indices = (:, :, grid.Nz),
                                                filename = near_global_filename(label, timestepper, filter, free_surface) * "_average",
                                                with_halos = true,
                                                overwrite_existing = true,
                                                array_type = Array{Float32})

    add_callback!(ocean, near_global_progress, progress_interval)

    cold_wall = @elapsed run!(ocean)

    ocean.Δt        = Δt
    ocean.stop_time = stop_time
    production_iter₀     = iteration(ocean)
    production_wall      = @elapsed run!(ocean)
    production_steps     = iteration(ocean) - production_iter₀

    wall_time  = cold_wall + production_wall
    iterations = iteration(ocean)
    seconds_per_step = production_steps == 0 ? NaN : production_wall / production_steps

    label = something(label, near_global_label(timestepper, filter, free_surface))
    @info @sprintf("[near_global] %-14s wall: %s (%d steps, cold %s) -> %.4f s/step",
                   label, prettytime(wall_time), iterations, prettytime(cold_wall), seconds_per_step)

    return (; ocean, label, wall_time, iterations, seconds_per_step)
end

function near_global_progress(sim)
    ocean = sim
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    step_time = 1e-9 * (time_ns() - TimestepperTestCases.wall_clock[])
    @info @sprintf("Iter: %d, time: %s, Δt: %s, max|u|: (%.2e, %.2e, %.2e), extrema(T): (%.2f, %.2f), wall: %s",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt),
                   maximum(abs, interior(u)), maximum(abs, interior(v)), maximum(abs, interior(w)),
                   maximum(interior(T)), minimum(interior(T)), prettytime(step_time))
    TimestepperTestCases.wall_clock[] = time_ns()
    return nothing
end

near_global_label(timestepper, filter, fs::SplitExplicitFreeSurface) = timestepper === :QuasiAdamsBashforth2 ? "AB2-SE-$(filter)" : "RK-SE-$(filter)"
near_global_label(timestepper, filter, fs::ImplicitFreeSurface) = timestepper === :QuasiAdamsBashforth2 ? "AB2-IM" : "RK-IM"

near_global_filename(label, timestepper, filter, fs) = "near_global_" * something(label, near_global_label(timestepper, filter, fs))

function run_near_global_cost(; arch = CPU(),
                                timestep = scheme -> near_global_timestep(Val(scheme)),
                                stop_time = 365days,
                                variants = near_global_discretizations())

    grid = near_global_grid(arch)
    results = []
    for d in variants
        result = near_global(d; arch, grid, Δt = timestep(d.timestepper), stop_time)
        push!(results, result)
    end

    @info "[near_global] cost summary"
    for r in results
        @info @sprintf("  %-14s  %8.4f s/step   (%s over %d steps)",
                       r.label, r.seconds_per_step, prettytime(r.wall_time), r.iterations)
    end

    return results
end

"""
    near_global(d::Discretization; arch = CPU(), kw...)

Run the near-global case with the discretization `d` of Table 1. Unlike the idealized cases the time step is
not derived from a theoretical limit here -- see `near_global_timestep` -- but the ratios between the schemes
are the same.
"""
function near_global(d::Discretization; arch = CPU(), grid = near_global_grid(arch), kw...)
    # `Discretization` carries the tracer advection of the idealized cases, which the near-global case overrides
    # with its own vertically implicit WENO7; only a scheme that departs from that shared default is forwarded, so
    # a discretization naming its own scheme reaches the model with it and every other one keeps the near-global.
    tracer_advection = d.tracer_advection === TimestepperTestCases.tracer_advection ? nothing : d.tracer_advection

    return near_global(d.timestepper; arch, grid,
                       free_surface = d.implicit_free_surface ? ImplicitFreeSurface() : nothing,
                       averaging_kernel = d.averaging_kernel,
                       barotropic_timestepper = d.barotropic_timestepper,
                       slow_forcing = d.slow_forcing,
                       tracer_advection,
                       label = d.label, kw...)
end
