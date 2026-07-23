using Oceananigans
using Oceananigans.Units
using NumericalEarth
using Dates: DateTime
using Printf

using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces:
    averaging_shape_function,
    LowDissipationAveragingKernel,
    SymmetricTrigAveragingKernel,
    WideTrig74AveragingKernel,
    WideTrig2AveragingKernel,
    ForwardBackwardScheme,
    RungeKutta3Scheme

using Oceananigans.BuoyancyFormulations: LinearEquationOfState

"""
    near_global_kernel(name::Symbol)

Return the split-explicit averaging kernel corresponding to a filter `name`.

Available kernels: `:SM05` (Shchepetkin--McWilliams shape function), `:mu2`
(`LowDissipationAveragingKernel`, ``\\mu_2 = 0``), `:trig` (`SymmetricTrigAveragingKernel`,
narrow window ``\\mu_2 = \\mu_3 = 0``), and the wide-window `:trig74` / `:trig2`
(`WideTrig74AveragingKernel` / `WideTrig2AveragingKernel`), the paper's default and
most stratification-robust filters.
"""
near_global_kernel(name::Symbol)   = near_global_kernel(Val(name))
near_global_kernel(::Val{:SM05})   = averaging_shape_function
near_global_kernel(::Val{:mu2})    = LowDissipationAveragingKernel()
near_global_kernel(::Val{:trig})   = SymmetricTrigAveragingKernel()
near_global_kernel(::Val{:trig74}) = WideTrig74AveragingKernel()
near_global_kernel(::Val{:trig2})  = WideTrig2AveragingKernel()

"""
    near_global_grid(arch; Nx, Ny, Nz, depth, latitude, longitude,
                     minimum_depth, interpolation_passes, major_basins)

Build the 1/4 degree near-global `ImmersedBoundaryGrid` with ETOPO bathymetry,
following the NumericalEarth `near_global_ocean_simulation` example.

$(SIGNATURES)
"""
function near_global_grid(arch = CPU();
                          Nx = 1440,
                          Ny = 600,
                          Nz = 70,
                          depth = 6000meters,
                          latitude  = (-75, 75),
                          longitude = (0, 360),
                          minimum_depth = 10meters,
                          interpolation_passes = 5,
                          major_basins = 3)

    # `mutable=true` gives a z-star (moving, free-surface-following) vertical coordinate,
    # consistent with the r/z-star formulation of the paper and the idealized cases.
    z = ExponentialDiscretization(Nz, -depth, 0, mutable=true)

    grid = LatitudeLongitudeGrid(arch;
                                 size = (Nx, Ny, Nz),
                                 halo = (7, 7, 7),
                                 z, latitude, longitude)

    bottom_height = regrid_bathymetry(grid; minimum_depth, interpolation_passes, major_basins)

    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)
end

"""
    near_global(timestepper::Symbol = :SplitRungeKutta3; kwargs...)

Set up and run the near-global quarter-degree test case, measuring wall-clock cost and
the numerical dissipation of the T and S variance.

$(SIGNATURES)

# Arguments
- `timestepper`: `:SplitRungeKutta3` (RK-SE / RK-IM) or `:QuasiAdamsBashforth2` (AB2-SE).

# Keyword Arguments
- `arch`: architecture, defaults to `CPU()`; pass `GPU()` (with CUDA loaded) for a real run.
- `free_surface`: pass an `ImplicitFreeSurface` for the RK-IM variant. When left as `nothing`,
                  a `SplitExplicitFreeSurface(grid; cfl, averaging_kernel)` is built (RK-SE / AB2-SE).
- `filter`: averaging-filter name for the split-explicit free surface (see `near_global_kernel`).
- `barotropic_timestepper`: barotropic sub-step scheme, `RungeKutta3Scheme()` (default) or `ForwardBackwardScheme()`.
- `cfl`: barotropic Courant number used to size the number of substeps (default `0.7`).
- `Δt`: baroclinic time step (default `25minutes`).
- `stop_time`: simulation length (default `365days`).
- `dissipation`: attach the `BuoyancyVarianceDissipation` diagnostic (default `true`).
- `equation_of_state`: defaults to a `LinearEquationOfState`, required by the buoyancy-dissipation diagnostic.
- `grid`: a prebuilt grid, to avoid re-deriving the bathymetry across variants.

# Returns
- A named tuple `(; simulation, ocean, label, wall_time, iterations, seconds_per_step)`.
"""
function near_global(timestepper::Symbol = :SplitRungeKutta3;
                     arch = CPU(),
                     grid = near_global_grid(arch),
                     free_surface = nothing,
                     filter::Symbol = :trig74,
                     barotropic_timestepper = RungeKutta3Scheme(),
                     cfl = 0.7,
                     Δt = 25minutes,
                     stop_time = 365days,
                     dissipation = true,
                     equation_of_state = LinearEquationOfState(),
                     label = nothing,
                     init_date = DateTime(1992, 1, 1),
                     progress_interval = TimeInterval(5days),
                     surface_output_interval = TimeInterval(1days),
                     dissipation_output_interval = TimeInterval(30days))

    if free_surface === nothing
        free_surface = SplitExplicitFreeSurface(grid; cfl,
                                                averaging_kernel = near_global_kernel(filter),
                                                timestepper = barotropic_timestepper)
    end

    ocean = ocean_simulation(grid; free_surface, timestepper, Δt, equation_of_state)

    set!(ocean.model, MetadataSet(:temperature, :salinity; dataset=ECCO4Monthly(), date=init_date))

    # Exact buoyancy variance-budget diagnostic (buoyancy is built from T and S through the
    # linear equation of state). NOTE: attach to the ocean simulation (not the coupled one) so
    # the callback receives `ocean.model`; verify the ocean callbacks fire under the coupled
    # `run!` on the target machine.
    if dissipation
        ϵb = BuoyancyVarianceDissipation(grid)
        add_callback!(ocean, ϵb, IterationInterval(1))

        diss = BuoyancyVarianceDissipationComputations.flatten_dissipation_fields(ϵb)
        ocean.output_writers[:dissipation] = JLD2Writer(ocean.model, diss;
                                                        schedule = dissipation_output_interval,
                                                        filename = _near_global_filename(label, timestepper, filter, free_surface) * "_dissipation",
                                                        overwrite_existing = true)
    end

    surface = merge(ocean.model.tracers, ocean.model.velocities)
    ocean.output_writers[:surface] = JLD2Writer(ocean.model, surface;
                                                schedule = surface_output_interval,
                                                indices = (:, :, grid.Nz),
                                                filename = _near_global_filename(label, timestepper, filter, free_surface) * "_surface",
                                                with_halos = true,
                                                overwrite_existing = true,
                                                array_type = Array{Float32})

    atmosphere = JRA55PrescribedAtmosphere(arch)
    radiation  = JRA55PrescribedRadiation(arch)
    land       = JRA55PrescribedLand(arch)
    coupled_model = OceanOnlyModel(ocean; atmosphere, land, radiation)

    simulation = Simulation(coupled_model; Δt, stop_time)
    add_callback!(simulation, near_global_progress, progress_interval)

    wall_time = @elapsed run!(simulation)
    iterations = iteration(simulation)
    seconds_per_step = iterations == 0 ? NaN : wall_time / iterations

    label = something(label, _near_global_label(timestepper, filter, free_surface))
    @info @sprintf("[near_global] %-14s wall time: %s over %d steps -> %.4f s/step",
                   label, prettytime(wall_time), iterations, seconds_per_step)

    return (; simulation, ocean, label, wall_time, iterations, seconds_per_step)
end

function near_global_progress(sim)
    ocean = sim.model.ocean
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

_free_surface_tag(::SplitExplicitFreeSurface) = "SE"
_free_surface_tag(::ImplicitFreeSurface) = "IM"

_near_global_label(timestepper, filter, fs::SplitExplicitFreeSurface) =
    timestepper === :QuasiAdamsBashforth2 ? "AB2-SE-$(filter)" : "RK-SE-$(filter)"
_near_global_label(timestepper, filter, fs::ImplicitFreeSurface) =
    timestepper === :QuasiAdamsBashforth2 ? "AB2-IM" : "RK-IM"

_near_global_filename(label, timestepper, filter, fs) =
    "near_global_" * something(label, _near_global_label(timestepper, filter, fs))

"""
    near_global_variants()

Return the list of `(; label, timestepper, filter, barotropic, implicit)` describing the
cost/dissipation comparison: AB2-SE (forward--backward barotropic), RK-SE with each averaging
filter (RK3 barotropic sub-step), and RK-IM.

$(SIGNATURES)
"""
function near_global_variants()
    RK3 = RungeKutta3Scheme()
    FB  = ForwardBackwardScheme()
    return [(; label = "AB2-SE",     timestepper = :QuasiAdamsBashforth2, filter = :SM05,   barotropic = FB,  implicit = false),
            (; label = "RK-SE-SM05", timestepper = :SplitRungeKutta3,     filter = :SM05,   barotropic = RK3, implicit = false),
            (; label = "RK-SE-mu2",  timestepper = :SplitRungeKutta3,     filter = :mu2,    barotropic = RK3, implicit = false),
            (; label = "RK-SE-trig74", timestepper = :SplitRungeKutta3,   filter = :trig74, barotropic = RK3, implicit = false),
            (; label = "RK-SE-trig2",  timestepper = :SplitRungeKutta3,   filter = :trig2,  barotropic = RK3, implicit = false),
            (; label = "RK-IM",      timestepper = :SplitRungeKutta3,     filter = :SM05,   barotropic = RK3, implicit = true)]
end

"""
    run_near_global_cost(; arch, Δt, Δt_ab2, stop_time, variants)

Run every variant of `near_global_variants` on a shared bathymetry grid and return a vector
of the per-variant cost results. Each variant also writes its T/S variance-dissipation output.

$(SIGNATURES)

The AB2 variant uses `Δt_ab2` (defaults to `Δt / 2`, reflecting its smaller stability range);
all Runge--Kutta variants use `Δt`.
"""
function run_near_global_cost(; arch = CPU(),
                                Δt = 25minutes,
                                Δt_ab2 = Δt / 2,
                                stop_time = 365days,
                                variants = near_global_variants())

    grid = near_global_grid(arch)
    results = []
    for v in variants
        free_surface = v.implicit ? ImplicitFreeSurface() : nothing
        Δt_v = v.timestepper === :QuasiAdamsBashforth2 ? Δt_ab2 : Δt
        result = near_global(v.timestepper; arch, grid, free_surface,
                             filter = v.filter, barotropic_timestepper = v.barotropic,
                             Δt = Δt_v, stop_time, label = v.label)
        push!(results, result)
    end

    @info "[near_global] cost summary"
    for r in results
        @info @sprintf("  %-14s  %8.4f s/step   (%s over %d steps)",
                       r.label, r.seconds_per_step, prettytime(r.wall_time), r.iterations)
    end

    return results
end
