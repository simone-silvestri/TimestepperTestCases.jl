using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znodes, φnode
using Oceananigans.Operators: intrinsic_vector, Δxᶜᶜᶜ, Δyᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: static_column_depthᶜᶜᵃ
using Oceananigans.OutputWriters: WindowedTimeAverage
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

global_ocean_kernel(name::Symbol)     = global_ocean_kernel(Val(name))
global_ocean_kernel(::Val{:SM05})     = averaging_shape_function
global_ocean_kernel(::Val{:mu2})      = LowDissipationAveragingKernel()
global_ocean_kernel(::Val{:trig})     = SymmetricTrigAveragingKernel()
global_ocean_kernel(::Val{:widetrig}) = WideTrigAveragingKernel()
global_ocean_kernel(::Val{:optasym})  = OptimizedAsymmetricAveragingKernel()

"""
    global_ocean_vertical_discretization(Nz = 100, depth = 5000meters, surface_spacing = 5meters)

Mutable `ExponentialDiscretization` of `Nz` levels over `depth` [m], with the e-folding scale solved for a top
cell of `surface_spacing` [m].
"""
function global_ocean_vertical_discretization(Nz = 100, depth = 5000meters, surface_spacing = 5meters)
    function top_spacing(scale)
        z = ExponentialDiscretization(Nz, -depth, 0; scale)
        return z[Nz+1] - z[Nz]
    end

    scale = find_zero(scale -> top_spacing(scale) - surface_spacing, (depth / Nz, 1000depth))

    return ExponentialDiscretization(Nz, -depth, 0; scale, mutable = true)
end


# measured, not derived: c₁ and Δx both vary over the domain, so the binding c₁ k needs the stratification field
global_ocean_timestep(::Val{:SplitRungeKutta3}) = 20minutes

# every other scheme scales from that reference by θ★(scheme)/θ★(WRK3), as in the idealized cases
global_ocean_timestep(::Val{scheme}) where scheme = baroclinic_timestep(scheme, global_ocean_timestep(Val(:SplitRungeKutta3)))

"""
    global_ocean_barotropic_rate(grid; g = Oceananigans.defaults.gravitational_acceleration)

Largest `√(gH) 2√(Δx⁻² + Δy⁻²)` over the columns of `grid` [s⁻¹], the rate that binds the barotropic sub-cycle.
"""
global_ocean_barotropic_rate(grid; g = Oceananigans.defaults.gravitational_acceleration) =
    maximum(compute!(Field(KernelFunctionOperation{Center, Center, Nothing}(column_barotropic_rate, grid, g))))

# `g` enters as an argument: the global `Oceananigans.defaults` lives in host memory and faults inside a GPU kernel
@inline column_barotropic_rate(i, j, k, grid, g) = barotropic_speed(static_column_depthᶜᶜᵃ(i, j, grid); g) *
                                                   staggered_wavenumber(Δxᶜᶜᶜ(i, j, k, grid), Δyᶜᶜᶜ(i, j, k, grid))

"""
    global_ocean_substeps(barotropic_scheme, grid, scheme; averaging_kernel, Δt)

Barotropic substep count that puts the substep Courant number `c₀ k Δτ` at the `barotropic_cfl` of
`barotropic_scheme`, with `c₀ k` the rate of `global_ocean_barotropic_rate`.
"""
function global_ocean_substeps(barotropic_scheme, grid, scheme = :SplitRungeKutta3;
                              averaging_kernel = OptimizedAsymmetricAveragingKernel(),
                              Δt = global_ocean_timestep(Val(scheme)))
    return barotropic_substeps(barotropic_scheme; Δt, averaging_kernel, rate = global_ocean_barotropic_rate(grid))
end

"""
    global_ocean_closure(FT = Oceananigans.defaults.FloatType)

The CATKE closure of this case: the vertically implicit `CATKEVerticalDiffusivity` with the bottom distance
coefficient of the shear mixing length set to `Cᵇ = 0.01`, against the 0.28 of `CATKEMixingLength`.
"""
global_ocean_closure(FT = Oceananigans.defaults.FloatType) =
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
    global_ocean_wind_stress(grid; latitudes, stresses, reference_density)

The eastward stress of `zonal_wind_stress` as the pair of surface momentum fluxes `(τˣ, τʸ)` [m² s⁻²],
projected onto the grid directions and negated, a top flux `J` contributing `-J/Δz` to the boundary cell.
"""
function global_ocean_wind_stress(grid; latitudes = default_wind_stress_latitudes,
                                       stresses = default_wind_stress_values,
                                       reference_density = 1026)
    τˣ = Field{Face, Center, Nothing}(grid)
    τʸ = Field{Center, Face, Nothing}(grid)

    launch!(architecture(grid), grid, :xy, _set_global_ocean_wind_stress!,
            τˣ, τʸ, grid, latitudes, stresses, reference_density)

    Oceananigans.BoundaryConditions.fill_halo_regions!(τˣ)
    Oceananigans.BoundaryConditions.fill_halo_regions!(τʸ)

    return τˣ, τʸ
end

# the rotation angle is cell-centered, so both components are evaluated at centers and written at their own
# staggered index: the half-cell offset is well below the scale on which this profile varies
@kernel function _set_global_ocean_wind_stress!(τˣ, τʸ, grid, latitudes, stresses, reference_density)
    i, j = @index(Global, NTuple)

    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    τ = - zonal_wind_stress(φ, latitudes, stresses) / reference_density

    τᵢ, τⱼ = intrinsic_vector(i, j, 1, grid, τ, zero(τ))

    @inbounds τˣ[i, j, 1] = τᵢ
    @inbounds τʸ[i, j, 1] = τⱼ
end

"""
    global_ocean_surface_restoring(grid, name; piston_velocity, dataset)

Surface restoring of `name`, `:temperature` or `:salinity`, towards the `dataset` climatology, as a top flux
through `SurfaceFluxRestoring` so that the interior buoyancy-variance budget stays untouched.

$(SIGNATURES)

`piston_velocity` is in m day⁻¹ and sets the rate through the thickness of the top cell. The WOA fields are
used as stored, in-situ temperature and practical salinity, the case running a `LinearEquationOfState`.
"""
function global_ocean_surface_restoring(grid, name; piston_velocity = 1/6, dataset = WOAMonthly())
    zF = Array(znodes(grid, Face()))
    surface_thickness = zF[end] - zF[end-1]
    rate = piston_velocity / (surface_thickness * days)

    metadata = Metadata(name; dataset)
    restoring = DatasetRestoring(metadata, architecture(grid); rate,
                                 time_indices_in_memory = length(metadata))

    return SurfaceFluxRestoring(restoring)
end

function global_ocean(timestepper::Symbol = :SplitRungeKutta3;
                     arch = CPU(),
                     grid = global_ocean_grid(arch),
                     free_surface = nothing,
                     filter::Symbol = :optasym,
                     averaging_kernel = global_ocean_kernel(filter),
                     barotropic_timestepper = ForwardBackwardScheme(),
                     slow_forcing = FrozenSlowForcing(),
                     tracer_advection = nothing,
                     momentum_advection = nothing,
                     Δt = global_ocean_timestep(Val(timestepper)),
                     nominal_timestep = Δt,
                     cold_start_Δt = Δt / 3,
                     cold_start_duration = 60days,
                     stop_time = 720days,
                     cost_stop_time = min(stop_time, 365days),
                     dissipation = true,
                     statistics = true,
                     closure = global_ocean_closure(),
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
                     prefix = "global_ocean",
                     output = true,
                     init_date = DateTime(1993, 1, 1),
                     progress_interval = TimeInterval(5days),
                     surface_output_interval = TimeInterval(1days),
                     dissipation_output_interval = AveragedTimeInterval(30days),
                     statistics_output_interval = AveragedTimeInterval(7days),
                     checkpoint_interval = 60days,
                     pickup = false,
                     dissipation_start_time = cost_stop_time)

    if free_surface === nothing
        # Fixed substep count, set from the substep integrator's own stability limit rather than from a cfl
        # target, so that the barotropic Courant number is a stated property of the run.
        substeps = global_ocean_substeps(barotropic_timestepper, grid; averaging_kernel, Δt = nominal_timestep)
        free_surface = SplitExplicitFreeSurface(grid; substeps, averaging_kernel,
                                                timestepper = barotropic_timestepper,
                                                slow_forcing)
    end
    
    # the campaign the finished cases belong to ran at cfl 1; a tighter threshold moves more vertical flux
    # into the implicit solve, which carries its own diffusivity and is the quantity being measured
    time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = parse(Float64, get(ENV, "M1_VERTICAL_CFL", "1")))

    momentum_advection = something(momentum_advection, split_momentum_advection(nothing, time_discretization))
    tracer_advection   = something(tracer_advection, weno_tracer_advection(9, time_discretization))

    wind_fluxes = if wind_stress
        τˣ, τʸ = global_ocean_wind_stress(grid; latitudes = wind_stress_latitudes,
                                        stresses = wind_stress_values, reference_density)
        (; u = FluxBoundaryCondition(τˣ), v = FluxBoundaryCondition(τʸ))
    else
        NamedTuple()
    end

    restoring_fluxes = if surface_restoring
        restoring(name) = global_ocean_surface_restoring(grid, name; piston_velocity = restoring_piston_velocity,
                                                        dataset = restoring_dataset)
        (; T = restoring(:temperature), S = restoring(:salinity))
    else
        NamedTuple()
    end

    additional_surface_fluxes = merge(wind_fluxes, restoring_fluxes)

    ocean = ocean_simulation(grid; free_surface, timestepper = campaign_timestepper(timestepper), Δt = cold_start_Δt, stop_time = cold_start_duration,
                             closure, equation_of_state, momentum_advection, tracer_advection,
                             bottom_drag_coefficient, bottom_drag_background_velocity,
                             additional_surface_fluxes)

    Tmetadata = Metadatum(:temperature, dataset=ECCO2Daily(), date=init_date)
    Smetadata = Metadatum(:salinity,    dataset=ECCO2Daily(), date=init_date)
    
    download_from_artifacts(metadata_path(Tmetadata))
    download_from_artifacts(metadata_path(Smetadata))

    set!(ocean.model, T=Tmetadata, S=Smetadata)

    parent(ocean.model.tracers.T) .= max.(parent(ocean.model.tracers.T), -1.8)

    filename = global_ocean_filename(prefix, label, timestepper, filter, free_surface)
    surface  = merge(ocean.model.tracers, ocean.model.velocities)

    if output
        ocean.output_writers[:surface] = JLD2Writer(ocean.model, surface;
                                                    schedule = surface_output_interval,
                                                    indices = (:, :, grid.Nz),
                                                    filename = filename * "_surface",
                                                    with_halos = true,
                                                    overwrite_files = true,
                                                    array_type = Array{Float32})

        ocean.output_writers[:average] = JLD2Writer(ocean.model, surface;
                                                    schedule = dissipation_output_interval,
                                                    indices = (:, :, grid.Nz),
                                                    filename = filename * "_average",
                                                    with_halos = true,
                                                    overwrite_files = true,
                                                    array_type = Array{Float32})
    end

    # Without a checkpointer a job that dies at hour ten of a fourteen hour run leaves nothing to resume from,
    # and this configuration runs long enough for that to be the common failure. `cleanup` keeps only the
    # newest, so the chain costs one checkpoint of disk rather than one per interval.
    if !isnothing(checkpoint_interval)
        ocean.output_writers[:checkpointer] = Checkpointer(ocean.model;
                                                           schedule = TimeInterval(checkpoint_interval),
                                                           prefix = filename * "_checkpoint",
                                                           overwrite_files = true,
                                                           cleanup = true)
    end

    add_callback!(ocean, global_ocean_progress, progress_interval)

    # A resumed run restores a state that is already past the cold start and past the segment the cost is
    # measured over, so both phases are skipped and their timings are not reported: `seconds_per_step` belongs
    # to the run that measured it, not to this one.
    if pickup
        ocean.Δt = Δt
        cold_wall, production_wall, production_steps = 0.0, 0.0, 0
    else
        cold_wall = @elapsed run!(ocean)

        ocean.Δt        = Δt
        ocean.stop_time = cost_stop_time
        first_production_iteration = iteration(ocean)
        production_wall  = @elapsed run!(ocean)
        production_steps = iteration(ocean) - first_production_iteration
    end

    # The budget's averaging windows are phased from the moment its writer is attached, so the join time is
    # what puts them on a given grid: joining at day 1095 with a 90-day window closes at 1185, 1275 and so on.
    if time(ocean) < dissipation_start_time
        ocean.stop_time = dissipation_start_time
        @info @sprintf("[global_ocean] %s running to the budget join at day %.0f", label,
                       dissipation_start_time / 86400)
        flush(stderr)
        run!(ocean; pickup)
    end

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

        diss = merge(TimestepperTestCases.BuoyancyVarianceDissipationComputations.flatten_dissipation_fields(ϵb), (; Gbx, Gby, Gbz))
        ocean.output_writers[:dissipation] = JLD2Writer(ocean.model, diss;
                                                        schedule = dissipation_output_interval,
                                                        filename = filename * "_dissipation",
                                                        overwrite_files = true)
    end

    # Float64 because the variance is recovered as ψ² - ψ̄², which cancels below Float32 precision in the abyss
    if output && statistics
        statistics_writer = JLD2Writer(ocean.model, global_ocean_statistics(ocean.model);
                                       schedule = statistics_output_interval,
                                       filename = filename * "_statistics",
                                       overwrite_files = true,
                                       array_type = Array{Float64})

        ocean.output_writers[:statistics] = initialize_schedules!(statistics_writer, ocean.model)
    end

    # `run!` takes one step even when the stop time is already reached
    ocean.stop_time = stop_time
    # `pickup = true` hands the newest checkpoint of this prefix to `run!`, which restores the model and the
    # clock before stepping; the writers above were built first so the budget is attached from the restart on.
    dissipation_wall = stop_time > cost_stop_time ? @elapsed(run!(ocean)) : 0.0

    wall_time  = cold_wall + production_wall + dissipation_wall
    iterations = iteration(ocean)
    seconds_per_step = production_steps == 0 ? NaN : production_wall / production_steps
    simulated_years_per_day = global_ocean_simulated_years_per_day(nominal_timestep, seconds_per_step)

    label = something(label, global_ocean_label(timestepper, filter, free_surface))
    @info @sprintf("[%s] %-14s wall: %s (%d steps, cold %s) -> %.4f s/step, %.3f SYPD at Δt = %s",
                   prefix, label, prettytime(wall_time), iterations, prettytime(cold_wall), seconds_per_step,
                   simulated_years_per_day, prettytime(nominal_timestep))

    substeps = free_surface isa SplitExplicitFreeSurface ? length(free_surface.substepping.averaging_weights) : 0

    save_global_ocean_cost(prefix, label; Δt, nominal_timestep, substeps, wall_time, cold_wall, production_wall, dissipation_wall,
                          iterations, production_steps, seconds_per_step, simulated_years_per_day,
                          cost_stop_time, stop_time)

    return (; ocean, label, wall_time, iterations, seconds_per_step, nominal_timestep, simulated_years_per_day)
end

"""
    global_ocean_statistics(model)

The fields whose time averages make up the `_statistics` output: the prognostic `u, v, T, S, e, η` and the
buoyancy `b`, their squares `u², v², T², S², e², η², b²`, from which the variance over an averaging window follows
as `ψ² - ψ̄²`, and the diffusive dissipation of buoyancy variance by CATKE `χᴷ = 2 κᶜ (∂z b)² V`, volume-integrated
like the advective production of the `_dissipation` output.
"""
function global_ocean_statistics(model)
    grid = model.grid
    u, v, _ = model.velocities
    T, S, e = model.tracers.T, model.tracers.S, model.tracers.e
    η = model.free_surface.displacement
    b = Oceananigans.Models.buoyancy_operation(model)

    means   = (; u, v, T, S, e, η, b)
    squares = NamedTuple{map(name -> Symbol(name, "²"), keys(means))}(map(ψ -> ψ^2, values(means)))

    V  = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Face), Oceananigans.Operators.volume, grid)
    κᶜ = model.closure_fields.κc
    χᴷ = 2 * κᶜ * ∂z(b)^2 * V

    return merge(means, squares, (; χᴷ))
end

"""
    initialize_schedules!(writer, model)

Start the schedule of `writer`, and the averaging windows of its time-averaged outputs, at the current time of
`model`. `run!` initializes schedules only at iteration 0, so a writer added to a simulation that has already run
otherwise counts its windows from t = 0.
"""
function initialize_schedules!(writer, model)
    Oceananigans.initialize!(writer.schedule, model)
    foreach(output -> initialize_output_schedule!(output, model), values(writer.outputs))
    return writer
end

initialize_output_schedule!(output, model) = nothing
initialize_output_schedule!(output::WindowedTimeAverage, model) = Oceananigans.initialize!(output.schedule, model)

"""
    save_global_ocean_cost(prefix, label; kw...)

Write the timings of a run to `<prefix>_<label>_cost.jld2`, beside its field output. Each keyword becomes one
entry of the file.
"""
function save_global_ocean_cost(prefix, label; kw...)
    filename = prefix * "_" * label * "_cost.jld2"

    JLD2.jldopen(filename, "w") do file
        file["label"] = label
        for (name, value) in pairs(kw)
            file[string(name)] = value
        end
    end

    return filename
end

free_surface_report(fs) = ""

function free_surface_report(fs::ImplicitFreeSurface)
    solver = fs.implicit_step_solver
    hasproperty(solver, :preconditioned_conjugate_gradient_solver) || return ""
    cg = solver.preconditioned_conjugate_gradient_solver
    preconditioner = isnothing(cg.preconditioner) ? "none" : nameof(typeof(cg.preconditioner))
    return @sprintf(", CG: %d/%d iters, |r|: %.3e, reltol %.1e, preconditioner %s",
                    cg.iteration, cg.maxiter, cg.residual_norm(cg.residual), cg.reltol, preconditioner)
end

function global_ocean_progress(sim)
    ocean = sim
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    step_time = 1e-9 * (time_ns() - TimestepperTestCases.wall_clock[])
    @info @sprintf("Iter: %d, time: %s, Δt: %s, max|u|: (%.2e, %.2e, %.2e), extrema(T): (%.2f, %.2f), wall: %s%s",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt),
                   maximum(abs, interior(u)), maximum(abs, interior(v)), maximum(abs, interior(w)),
                   maximum(interior(T)), minimum(interior(T)), prettytime(step_time), free_surface_report(ocean.model.free_surface))

    # Julia block-buffers stderr when it is a file, so a batch log only appears when the process exits
    flush(stderr)

    TimestepperTestCases.wall_clock[] = time_ns()
    return nothing
end

global_ocean_label(timestepper, filter, fs::SplitExplicitFreeSurface) = timestepper === :QuasiAdamsBashforth2 ? "AB2-SE-$(filter)" : "RK-SE-$(filter)"
global_ocean_label(timestepper, filter, fs::ImplicitFreeSurface) = timestepper === :QuasiAdamsBashforth2 ? "AB2-IM" : "RK-IM"

global_ocean_filename(prefix, label, timestepper, filter, fs) = prefix * "_" * something(label, global_ocean_label(timestepper, filter, fs))


"""
    global_ocean(d::Discretization; arch = CPU(), kw...)

Run the global case with the discretization `d` of Table 1.
"""
function global_ocean(d::Discretization; arch = CPU(), grid = global_ocean_grid(arch), kw...)
    # only a tracer advection that departs from the shared default is forwarded; every other discretization
    # keeps the vertically implicit WENO9 of this case
    return global_ocean(d.timestepper; arch, grid,
                       free_surface = d.implicit_free_surface ? ImplicitFreeSurface() : nothing,
                       averaging_kernel = d.averaging_kernel,
                       barotropic_timestepper = d.barotropic_timestepper,
                       slow_forcing = d.slow_forcing,
                       tracer_advection = forwarded_tracer_advection(d),
                       label = d.label, kw...)
end
