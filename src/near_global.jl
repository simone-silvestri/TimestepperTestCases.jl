using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing, znodes
using NumericalEarth
using Dates: DateTime
using Printf

using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces:
    averaging_shape_function,
    LowDissipationAveragingKernel,
    SymmetricTrigAveragingKernel,
    WideTrigAveragingKernel,
    OptimizedAsymmetricAveragingKernel,
    ForwardBackwardScheme,
    RungeKutta3Scheme

using Oceananigans.BuoyancyFormulations: LinearEquationOfState
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation

using Downloads: Downloads
using JLD2
using Oceananigans.Architectures: architecture
using NumericalEarth.DataWrangling: metadata_path, DatasetRestoring, SurfaceFluxRestoring
using NumericalEarth.DataWrangling.WOA: WOAMonthly

# Loads NumericalEarthWOAExt, which carries the download method for the World Ocean Atlas files
using WorldOceanAtlasTools

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

near_global_kernel(name::Symbol)     = near_global_kernel(Val(name))
near_global_kernel(::Val{:SM05})     = averaging_shape_function
near_global_kernel(::Val{:mu2})      = LowDissipationAveragingKernel()
near_global_kernel(::Val{:trig})     = SymmetricTrigAveragingKernel()
near_global_kernel(::Val{:widetrig}) = WideTrigAveragingKernel()
near_global_kernel(::Val{:optasym})  = OptimizedAsymmetricAveragingKernel()

function near_global_grid(arch = CPU();
                          Nx = 1440,
                          Ny = 600,
                          Nz = 50,
                          depth = 5000meters,
                          latitude  = (-75, 75),
                          longitude = (0, 360),
                          # A 15 m floor leaves one-cell pinnacles standing among their 1800 m neighbours --
                          # the Solomon Sea has a 16.7 m column beside a dry cell -- and the flow squeezing
                          # through the constriction accelerates without bound. 30 m puts the floor two cell
                          # interfaces deeper, and the extra interpolation passes coarsen the native
                          # bathymetry gradually enough that an isolated spike is averaged away instead of
                          # aliased onto a single column.
                          minimum_depth = 30meters,
                          interpolation_passes = 20,
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

#####
##### Idealized wind stress, the profile of `WenoNeverworld`
#####

# Latitudes and zonal-mean zonal stresses in N m⁻², the `default_φs` and `default_τs` of WenoNeverworld:
# Southern Ocean westerlies at 0.2, trade easterlies at -0.1 in both hemispheres, a weak equatorial -0.02, and
# northern westerlies at 0.1 -- half the Southern Ocean, which is the asymmetry of the observed zonal mean.
const default_wind_stress_latitudes = (-70, -45, -15, 0, 15, 45, 70)
const default_wind_stress_values    = (0.0, 0.2, -0.1, -0.02, -0.1, 0.1, 0.0)

"""
    zonal_wind_stress(φ, latitudes, stresses)

Zonal-mean zonal wind stress at latitude `φ`, in N m⁻², interpolated between the knots `(latitudes, stresses)`.

Piecewise cubic between adjacent knots with vanishing derivative at each of them, as in `WenoNeverworld`, so
that the profile is smooth and every knot is a local extremum of the stress. Contrarily to `WenoNeverworld`
the domain reaches beyond the outermost knots, where the stress is held at zero -- continuous, the first and
last knot carrying no stress.
"""
@inline function zonal_wind_stress(φ, latitudes = default_wind_stress_latitudes,
                                      stresses = default_wind_stress_values)
    φ ≤ first(latitudes) && return zero(φ)
    φ ≥ last(latitudes)  && return zero(φ)

    k = findfirst(≥(φ), latitudes)
    φ₁, φ₂ = latitudes[k-1], latitudes[k]
    τ₁, τ₂ = stresses[k-1], stresses[k]

    t = (φ - φ₁) / (φ₂ - φ₁)

    return τ₁ * (2t^3 - 3t^2 + 1) + τ₂ * (3t^2 - 2t^3)
end

"""
    near_global_wind_stress(grid; latitudes, stresses, reference_density)

The zonal momentum flux of [`zonal_wind_stress`](@ref) as a surface field, ready to be passed to
`ocean_simulation` through `additional_surface_fluxes`.

Ocean-only, the surface momentum flux that `ocean_simulation` allocates is never filled -- it is the coupled
interface that writes it -- so without this the configuration has no momentum sink at the surface at all, and
the geostrophic adjustment from the initial state has nothing acting against it.

The stress is negated: a top flux `J` contributes `-J/Δz` to the tendency of the boundary cell, so an eastward
stress reaches the model as a negative flux. `WenoNeverworld` carries the same sign for the same reason.
"""
function near_global_wind_stress(grid; latitudes = default_wind_stress_latitudes,
                                       stresses = default_wind_stress_values,
                                       reference_density = 1026)
    τˣ = Field{Face, Center, Nothing}(grid)
    set!(τˣ, (λ, φ) -> - zonal_wind_stress(φ, latitudes, stresses) / reference_density)
    Oceananigans.BoundaryConditions.fill_halo_regions!(τˣ)

    return τˣ
end

"""
    near_global_surface_restoring(grid, name; piston_velocity, dataset)

Surface restoring of `name` -- `:temperature` or `:salinity` -- towards the `dataset` climatology, following
the `salinity_surface_restoring` of the OMIP configuration.

$(SIGNATURES)

The target is the monthly World Ocean Atlas, so the restoring carries the seasonal cycle: `WOAMonthly` is a
twelve-month climatology cycled by the `Cyclical` time indexing of `DatasetRestoring`, which puts the surface
back where the season says it should be rather than holding it at one month all year. It is also a smooth
climatology rather than a state estimate, so it constrains the large-scale surface without pulling the model's
eddy field towards a different realization of the eddies.

`SurfaceFluxRestoring` evaluates the restoring at the top cell alone and converts the tendency into a top
flux, `-G Δz`. That matters for this case: a restoring term is a buoyancy source that the buoyancy-variance
budget does not account for, so a restoring reaching the interior would bias the numerical diffusivity the case
exists to measure. Entering as a surface flux it leaves the interior budget untouched.

`piston_velocity` is in m day⁻¹ and sets the rate through the thickness of the top cell, so the restoring
timescale follows the vertical grid instead of being quoted independently of it. Temperature is restored as
well as salinity, contrarily to OMIP, because this configuration carries no atmosphere: with no bulk heat flux
the restoring is the only surface constraint on the temperature.

The WOA fields are used as they are stored -- in-situ temperature and practical salinity -- with none of the
TEOS-10 conversion the OMIP configuration applies, this case running the `LinearEquationOfState` in which the
two enter through constant expansion coefficients.
"""
function near_global_surface_restoring(grid, name; piston_velocity = 1/6, dataset = WOAMonthly())
    zF = Array(znodes(grid, Face()))
    surface_thickness = zF[end] - zF[end-1]
    rate = piston_velocity / (surface_thickness * days)

    metadata = Metadata(name; dataset)
    restoring = DatasetRestoring(metadata, architecture(grid); rate,
                                 time_indices_in_memory = length(metadata))

    return SurfaceFluxRestoring(restoring)
end

function near_global(timestepper::Symbol = :SplitRungeKutta3;
                     arch = CPU(),
                     grid = near_global_grid(arch),
                     free_surface = nothing,
                     filter::Symbol = :optasym,
                     averaging_kernel = near_global_kernel(filter),
                     barotropic_timestepper = ForwardBackwardScheme(),
                     slow_forcing = FrozenSlowForcing(),
                     tracer_advection = nothing,
                     momentum_advection = nothing,
                     boundary_scheme = nothing,
                     tracer_boundary_scheme = something(boundary_scheme, default_tracer_boundary_scheme),
                     momentum_boundary_scheme = something(boundary_scheme, default_momentum_boundary_scheme),
                     # ∇ref Δ ≈ 0.14 K and ≈0.033 g/kg on the quarter-degree grid, about a decade below the
                     # typical horizontal gradient, so a topographic step fires the constant candidate while
                     # smooth data keeps third order. The vertical scales are read off the stencil.
                     horizontal_temperature_reference_gradient = 5e-6,
                     vertical_temperature_reference_gradient = 0,
                     horizontal_salinity_reference_gradient = 1.2e-6,
                     vertical_salinity_reference_gradient = 0,
                     vertical_momentum_reference_gradient = 0,
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
                     reference_density = 1026,
                     wind_stress = true,
                     wind_stress_latitudes = default_wind_stress_latitudes,
                     wind_stress_values = default_wind_stress_values,
                     surface_restoring = true,
                     restoring_piston_velocity = 1/6,
                     restoring_dataset = WOAMonthly(),
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

    tracer_boundary_scheme   = boundary_scheme_value(tracer_boundary_scheme)
    momentum_boundary_scheme = boundary_scheme_value(momentum_boundary_scheme)

    # The horizontal momentum terms reconstruct a vorticity, a divergence flux and a squared velocity, so no
    # single reference gradient carries their units: there the oscillation scale is read off the stencil.
    momentum_advection = something(momentum_advection,
                                   split_momentum_advection(nothing, time_discretization, momentum_boundary_scheme,
                                                            0, vertical_momentum_reference_gradient))

    # Temperature and salinity carry their own reference gradients, so each takes its own scheme.
    vertically_implicit_weno7 = WENO(order=7; time_discretization)

    tracer_advection = something(tracer_advection,
                                 (T = tracer_advection_scheme(tracer_boundary_scheme, vertically_implicit_weno7;
                                                              time_discretization,
                                                              horizontal_reference_gradient = horizontal_temperature_reference_gradient,
                                                              vertical_reference_gradient = vertical_temperature_reference_gradient),
                                  S = tracer_advection_scheme(tracer_boundary_scheme, vertically_implicit_weno7;
                                                              time_discretization,
                                                              horizontal_reference_gradient = horizontal_salinity_reference_gradient,
                                                              vertical_reference_gradient = vertical_salinity_reference_gradient)))

    # A boundary scheme that departs from the default names the run it produces, so that its output does not
    # overwrite the default one it is meant to be compared against. The label is materialized only where there
    # is a suffix to append, `nothing` reaching `near_global_filename` as the request for the derived name.
    boundary_suffix = boundary_scheme_suffix(tracer_boundary_scheme, momentum_boundary_scheme)

    if !isempty(boundary_suffix)
        label = something(label, near_global_label(timestepper, filter, free_surface)) * boundary_suffix
    end

    wind_fluxes = if wind_stress
        τˣ = near_global_wind_stress(grid; latitudes = wind_stress_latitudes,
                                     stresses = wind_stress_values, reference_density)
        (; u = FluxBoundaryCondition(τˣ))
    else
        NamedTuple()
    end

    restoring_fluxes = if surface_restoring
        restoring(name) = near_global_surface_restoring(grid, name; piston_velocity = restoring_piston_velocity,
                                                        dataset = restoring_dataset)
        (; T = restoring(:temperature), S = restoring(:salinity))
    else
        NamedTuple()
    end

    additional_surface_fluxes = merge(wind_fluxes, restoring_fluxes)

    ocean = ocean_simulation(grid; free_surface, timestepper, Δt = cold_start_Δt, stop_time = cold_start_duration,
                             closure, equation_of_state, momentum_advection, tracer_advection,
                             bottom_drag_coefficient, bottom_drag_background_velocity,
                             additional_surface_fluxes)

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

    substeps = free_surface isa SplitExplicitFreeSurface ? length(free_surface.substepping.averaging_weights) : 0

    save_near_global_cost(label; Δt, substeps, wall_time, cold_wall, production_wall,
                          iterations, production_steps, seconds_per_step, stop_time)

    return (; ocean, label, wall_time, iterations, seconds_per_step)
end

"""
    save_near_global_cost(label; kw...)

Write the timings of a near-global run to `near_global_<label>_cost.jld2`, beside its field output, so that
the cost comparison of Section 5.4.1 can be assembled from the output directory rather than from the return
value of the runner. Each keyword becomes one entry of the file.
"""
function save_near_global_cost(label; kw...)
    filename = "near_global_" * label * "_cost.jld2"

    JLD2.jldopen(filename, "w") do file
        file["label"] = label
        for (name, value) in pairs(kw)
            file[string(name)] = value
        end
    end

    return filename
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

    # Julia block-buffers stderr when it is a file, so a batch log only appears when the process exits
    flush(stderr)

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

    results = []

    # A fresh grid per case: the mutable vertical coordinate keeps ηⁿ and the σ scalings on the grid, so a
    # shared one carries the surface state -- a NaN included -- from each case into the next
    for d in variants
        grid = near_global_grid(arch)
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
    return near_global(d.timestepper; arch, grid,
                       free_surface = d.implicit_free_surface ? ImplicitFreeSurface() : nothing,
                       averaging_kernel = d.averaging_kernel,
                       barotropic_timestepper = d.barotropic_timestepper,
                       slow_forcing = d.slow_forcing,
                       tracer_advection = forwarded_tracer_advection(d),
                       label = d.label, kw...)
end
