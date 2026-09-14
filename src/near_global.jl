using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znodes, φnode
using Oceananigans.Operators: intrinsic_vector, Δxᶜᶜᶜ, Δyᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: static_column_depthᶜᶜᵃ
using Oceananigans.Utils: launch!
using Roots: find_zero
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

"""
    near_global_vertical_discretization(Nz = 100, depth = 5000meters, surface_spacing = 5meters)

Mutable `ExponentialDiscretization` of `Nz` levels over `depth` [m], with the e-folding scale solved for a top
cell of `surface_spacing` [m].
"""
function near_global_vertical_discretization(Nz = 100, depth = 5000meters, surface_spacing = 5meters)
    function top_spacing(scale)
        z = ExponentialDiscretization(Nz, -depth, 0; scale)
        return z[Nz+1] - z[Nz]
    end

    scale = find_zero(scale -> top_spacing(scale) - surface_spacing, (depth / Nz, 1000depth))

    return ExponentialDiscretization(Nz, -depth, 0; scale, mutable = true)
end

function near_global_grid(arch = CPU();
                          Nx = 1440,
                          Ny = 600,
                          Nz = 50,
                          depth = 5000meters,
                          latitude  = (-60, 60),
                          longitude = (0, 360),
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

# measured, not derived: c₁ and Δx both vary over the domain, so the binding c₁ k needs the stratification field
near_global_timestep(::Val{:SplitRungeKutta3}) = 20minutes

# every other scheme scales from that reference by θ★(scheme)/θ★(WRK3), as in the idealized cases
near_global_timestep(::Val{scheme}) where scheme = baroclinic_timestep(scheme, near_global_timestep(Val(:SplitRungeKutta3)))

"""
    near_global_barotropic_rate(grid; g = Oceananigans.defaults.gravitational_acceleration)

Largest `√(gH) 2√(Δx⁻² + Δy⁻²)` over the columns of `grid` [s⁻¹], the rate that binds the barotropic sub-cycle.
"""
near_global_barotropic_rate(grid; g = Oceananigans.defaults.gravitational_acceleration) =
    maximum(compute!(Field(KernelFunctionOperation{Center, Center, Nothing}(column_barotropic_rate, grid, g))))

# `g` enters as an argument: the global `Oceananigans.defaults` lives in host memory and faults inside a GPU kernel
@inline column_barotropic_rate(i, j, k, grid, g) = barotropic_speed(static_column_depthᶜᶜᵃ(i, j, grid); g) *
                                                   staggered_wavenumber(Δxᶜᶜᶜ(i, j, k, grid), Δyᶜᶜᶜ(i, j, k, grid))

"""
    near_global_substeps(barotropic_scheme, grid, scheme; averaging_kernel, Δt)

Barotropic substep count that puts the substep Courant number `c₀ k Δτ` at the `barotropic_cfl` of
`barotropic_scheme`, with `c₀ k` the rate of `near_global_barotropic_rate`.
"""
function near_global_substeps(barotropic_scheme, grid, scheme = :SplitRungeKutta3;
                              averaging_kernel = OptimizedAsymmetricAveragingKernel(),
                              Δt = near_global_timestep(Val(scheme)))
    return barotropic_substeps(barotropic_scheme; Δt, averaging_kernel, rate = near_global_barotropic_rate(grid))
end

"""
    near_global_discretizations()

The cases of Table 1 run on the near-global configuration: `discretizations()` without `WRK3-UP`.
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

Zonal-mean zonal wind stress at latitude `φ` [N m⁻²], piecewise cubic between the knots `(latitudes, stresses)`
with vanishing derivative at each of them, and zero beyond the outermost knots.
"""
@inline function zonal_wind_stress(φ, latitudes = default_wind_stress_latitudes,
                                      stresses = default_wind_stress_values)
    φ ≤ first(latitudes) && return zero(φ)
    φ ≥ last(latitudes)  && return zero(φ)

    # a bounded loop, since the `nothing` branch of `findfirst` does not compile in a GPU kernel
    k = 2
    while k < length(latitudes) && latitudes[k] < φ
        k += 1
    end

    φ₁, φ₂ = latitudes[k-1], latitudes[k]
    τ₁, τ₂ = stresses[k-1], stresses[k]

    t = (φ - φ₁) / (φ₂ - φ₁)

    return τ₁ * (2t^3 - 3t^2 + 1) + τ₂ * (3t^2 - 2t^3)
end

"""
    near_global_wind_stress(grid; latitudes, stresses, reference_density)

The eastward stress of `zonal_wind_stress` as the pair of surface momentum fluxes `(τˣ, τʸ)` [m² s⁻²],
projected onto the grid directions and negated, a top flux `J` contributing `-J/Δz` to the boundary cell.
"""
function near_global_wind_stress(grid; latitudes = default_wind_stress_latitudes,
                                       stresses = default_wind_stress_values,
                                       reference_density = 1026)
    τˣ = Field{Face, Center, Nothing}(grid)
    τʸ = Field{Center, Face, Nothing}(grid)

    launch!(architecture(grid), grid, :xy, _set_near_global_wind_stress!,
            τˣ, τʸ, grid, latitudes, stresses, reference_density)

    Oceananigans.BoundaryConditions.fill_halo_regions!(τˣ)
    Oceananigans.BoundaryConditions.fill_halo_regions!(τʸ)

    return τˣ, τʸ
end

# the rotation angle is cell-centered, so both components are evaluated at centers and written at their own
# staggered index: the half-cell offset is well below the scale on which this profile varies
@kernel function _set_near_global_wind_stress!(τˣ, τʸ, grid, latitudes, stresses, reference_density)
    i, j = @index(Global, NTuple)

    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    τ = - zonal_wind_stress(φ, latitudes, stresses) / reference_density

    τᵢ, τⱼ = intrinsic_vector(i, j, 1, grid, τ, zero(τ))

    @inbounds τˣ[i, j, 1] = τᵢ
    @inbounds τʸ[i, j, 1] = τⱼ
end

"""
    near_global_surface_restoring(grid, name; piston_velocity, dataset)

Surface restoring of `name`, `:temperature` or `:salinity`, towards the `dataset` climatology, as a top flux
through `SurfaceFluxRestoring` so that the interior buoyancy-variance budget stays untouched.

$(SIGNATURES)

`piston_velocity` is in m day⁻¹ and sets the rate through the thickness of the top cell. The WOA fields are
used as stored, in-situ temperature and practical salinity, the case running a `LinearEquationOfState`.
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
                     Δt = near_global_timestep(Val(timestepper)),
                     nominal_timestep = Δt,
                     cold_start_Δt = Δt / 3,
                     cold_start_duration = 60days,
                     stop_time = 720days,
                     cost_stop_time = min(stop_time, 365days),
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
                     prefix = "near_global",
                     output = true,
                     init_date = DateTime(1993, 1, 1),
                     progress_interval = TimeInterval(5days),
                     surface_output_interval = TimeInterval(1days),
                     dissipation_output_interval = AveragedTimeInterval(30days))

    if free_surface === nothing
        # Fixed substep count, set from the substep integrator's own stability limit rather than from a cfl
        # target, so that the barotropic Courant number is a stated property of the run.
        substeps = near_global_substeps(barotropic_timestepper, grid; averaging_kernel, Δt = nominal_timestep)
        free_surface = SplitExplicitFreeSurface(grid; substeps, averaging_kernel,
                                                timestepper = barotropic_timestepper,
                                                slow_forcing)
    end
    
    time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5)

    tracer_boundary_scheme   = boundary_scheme_value(tracer_boundary_scheme)
    momentum_boundary_scheme = boundary_scheme_value(momentum_boundary_scheme)

    momentum_advection = something(momentum_advection,
                                   split_momentum_advection(nothing, time_discretization, momentum_boundary_scheme))

    vertically_implicit_weno7 = WENO(order=7; time_discretization)

    tracer_advection = something(tracer_advection,
                                 tracer_advection_scheme(tracer_boundary_scheme, vertically_implicit_weno7;
                                                         time_discretization))

    # a boundary scheme that departs from the default names the run it produces, so that its output does not
    # overwrite the default one
    boundary_suffix = boundary_scheme_suffix(tracer_boundary_scheme, momentum_boundary_scheme)

    if !isempty(boundary_suffix)
        label = something(label, near_global_label(timestepper, filter, free_surface)) * boundary_suffix
    end

    wind_fluxes = if wind_stress
        τˣ, τʸ = near_global_wind_stress(grid; latitudes = wind_stress_latitudes,
                                        stresses = wind_stress_values, reference_density)
        (; u = FluxBoundaryCondition(τˣ), v = FluxBoundaryCondition(τʸ))
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

    filename = near_global_filename(prefix, label, timestepper, filter, free_surface)
    surface  = merge(ocean.model.tracers, ocean.model.velocities)

    if output
        ocean.output_writers[:surface] = JLD2Writer(ocean.model, surface;
                                                    schedule = surface_output_interval,
                                                    indices = (:, :, grid.Nz),
                                                    filename = filename * "_surface",
                                                    with_halos = true,
                                                    overwrite_existing = true,
                                                    array_type = Array{Float32})

        ocean.output_writers[:average] = JLD2Writer(ocean.model, surface;
                                                    schedule = dissipation_output_interval,
                                                    indices = (:, :, grid.Nz),
                                                    filename = filename * "_average",
                                                    with_halos = true,
                                                    overwrite_existing = true,
                                                    array_type = Array{Float32})
    end

    add_callback!(ocean, near_global_progress, progress_interval)

    cold_wall = @elapsed run!(ocean)

    ocean.Δt        = Δt
    ocean.stop_time = cost_stop_time
    first_production_iteration = iteration(ocean)
    production_wall  = @elapsed run!(ocean)
    production_steps = iteration(ocean) - first_production_iteration

    # the buoyancy-variance budget runs every iteration, so it joins only after the segment `seconds_per_step`
    # is measured over
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
                                                        filename = filename * "_dissipation",
                                                        overwrite_existing = true)
    end

    # `run!` takes one step even when the stop time is already reached
    ocean.stop_time = stop_time
    dissipation_wall = stop_time > cost_stop_time ? @elapsed(run!(ocean)) : 0.0

    wall_time  = cold_wall + production_wall + dissipation_wall
    iterations = iteration(ocean)
    seconds_per_step = production_steps == 0 ? NaN : production_wall / production_steps
    simulated_years_per_day = near_global_simulated_years_per_day(nominal_timestep, seconds_per_step)

    label = something(label, near_global_label(timestepper, filter, free_surface))
    @info @sprintf("[%s] %-14s wall: %s (%d steps, cold %s) -> %.4f s/step, %.3f SYPD at Δt = %s",
                   prefix, label, prettytime(wall_time), iterations, prettytime(cold_wall), seconds_per_step,
                   simulated_years_per_day, prettytime(nominal_timestep))

    substeps = free_surface isa SplitExplicitFreeSurface ? length(free_surface.substepping.averaging_weights) : 0

    save_near_global_cost(prefix, label; Δt, nominal_timestep, substeps, wall_time, cold_wall, production_wall, dissipation_wall,
                          iterations, production_steps, seconds_per_step, simulated_years_per_day,
                          cost_stop_time, stop_time)

    return (; ocean, label, wall_time, iterations, seconds_per_step, nominal_timestep, simulated_years_per_day)
end

"""
    near_global_simulated_years_per_day(nominal_timestep, seconds_per_step)

Simulated years per wall-clock day of a run at its production time step `nominal_timestep` [s], whatever step
`seconds_per_step` was measured at.
"""
near_global_simulated_years_per_day(nominal_timestep, seconds_per_step) = nominal_timestep / seconds_per_step / 365

"""
    save_near_global_cost(prefix, label; kw...)

Write the timings of a run to `<prefix>_<label>_cost.jld2`, beside its field output. Each keyword becomes one
entry of the file.
"""
function save_near_global_cost(prefix, label; kw...)
    filename = prefix * "_" * label * "_cost.jld2"

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

near_global_filename(prefix, label, timestepper, filter, fs) = prefix * "_" * something(label, near_global_label(timestepper, filter, fs))

function run_near_global_cost(; arch = CPU(),
                                case = near_global,
                                timestep = scheme -> near_global_timestep(Val(scheme)),
                                stop_time = 365days,
                                variants = near_global_discretizations())

    results = []

    # each case builds its own grid: the mutable vertical coordinate keeps ηⁿ and the σ scalings on the grid,
    # so a shared one carries the surface state from one case into the next
    for d in variants
        push!(results, case(d; arch, Δt = timestep(d.timestepper), stop_time, dissipation = false))
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

Run the near-global case with the discretization `d` of Table 1.
"""
function near_global(d::Discretization; arch = CPU(), grid = near_global_grid(arch), kw...)
    # only a tracer advection that departs from the shared default is forwarded; every other discretization
    # keeps the vertically implicit WENO7 of this case
    return near_global(d.timestepper; arch, grid,
                       free_surface = d.implicit_free_surface ? ImplicitFreeSurface() : nothing,
                       averaging_kernel = d.averaging_kernel,
                       barotropic_timestepper = d.barotropic_timestepper,
                       slow_forcing = d.slow_forcing,
                       tracer_advection = forwarded_tracer_advection(d),
                       label = d.label, kw...)
end
