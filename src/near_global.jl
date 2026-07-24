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

near_global_kernel(name::Symbol)   = near_global_kernel(Val(name))
near_global_kernel(::Val{:SM05})   = averaging_shape_function
near_global_kernel(::Val{:mu2})    = LowDissipationAveragingKernel()
near_global_kernel(::Val{:trig})   = SymmetricTrigAveragingKernel()
near_global_kernel(::Val{:trig74}) = WideTrig74AveragingKernel()
near_global_kernel(::Val{:trig2})  = WideTrig2AveragingKernel()

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

near_global_timestep(::Val{:QuasiAdamsBashforth2}) = 8.5minutes
near_global_timestep(::Val{:SplitRungeKutta3}) = 25minutes

function near_global(timestepper::Symbol = :SplitRungeKutta3;
                     arch = CPU(),
                     grid = near_global_grid(arch),
                     free_surface = nothing,
                     filter::Symbol = :trig74,
                     barotropic_timestepper = ForwardBackwardScheme(),
                     cfl = 0.7,
                     Δt = near_global_timestep(Val(timestepper)),
                     cold_start_Δt = Δt / 3,
                     cold_start_duration = 60days,
                     stop_time = 720days,
                     dissipation = true,
                     equation_of_state = LinearEquationOfState(),
                     label = nothing,
                     init_date = DateTime(1993, 1, 1),
                     progress_interval = TimeInterval(5days),
                     surface_output_interval = TimeInterval(1days),
                     dissipation_output_interval = AveragedTimeInterval(30days))

    if free_surface === nothing
        free_surface = SplitExplicitFreeSurface(grid; cfl,
                                                fixed_Δt = Δt + 2minutes,
                                                averaging_kernel = near_global_kernel(filter),
                                                timestepper = barotropic_timestepper)
    end
    
    time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5)
    momentum_advection = WENOVectorInvariant(; time_discretization)
    tracer_advection = WENO(order=7; minimum_buffer_upwind_order=3, time_discretization)

    ocean = ocean_simulation(grid; free_surface, timestepper, Δt, equation_of_state, momentum_advection, tracer_advection)

    Tmetadata = Metadatum(:temperature, dataset=ECCO2Daily(), date=init_date)
    Smetadata = Metadatum(:salinity,    dataset=ECCO2Daily(), date=init_date)
    
    download_from_artifacts(metadata_path(Tmetadata))
    download_from_artifacts(metadata_path(Smetadata))

    set!(ocean.model, T=Tmetadata, S=Smetadata)

    if dissipation
        ϵb = BuoyancyVarianceDissipation(grid)
        add_callback!(ocean, ϵb, IterationInterval(1))

        diss = BuoyancyVarianceDissipationComputations.flatten_dissipation_fields(ϵb)
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
                                                filename = near_global_filename(label, timestepper, filter, free_surface) * "_surface",
                                                with_halos = true,
                                                overwrite_existing = true,
                                                array_type = Array{Float32})

    atmosphere = JRA55PrescribedAtmosphere(arch)
    radiation  = JRA55PrescribedRadiation(arch)
    land       = JRA55PrescribedLand(arch)
    coupled_model = OceanOnlyModel(ocean; atmosphere, land, radiation)

    simulation = Simulation(coupled_model; Δt = cold_start_Δt, stop_time = cold_start_duration)
    add_callback!(simulation, near_global_progress, progress_interval)

    cold_wall = @elapsed run!(simulation)

    simulation.Δt        = Δt
    simulation.stop_time = stop_time
    production_iter₀     = iteration(simulation)
    production_wall      = @elapsed run!(simulation)
    production_steps     = iteration(simulation) - production_iter₀

    wall_time  = cold_wall + production_wall
    iterations = iteration(simulation)
    seconds_per_step = production_steps == 0 ? NaN : production_wall / production_steps

    label = something(label, near_global_label(timestepper, filter, free_surface))
    @info @sprintf("[near_global] %-14s wall: %s (%d steps, cold %s) -> %.4f s/step",
                   label, prettytime(wall_time), iterations, prettytime(cold_wall), seconds_per_step)

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

near_global_label(timestepper, filter, fs::SplitExplicitFreeSurface) = timestepper === :QuasiAdamsBashforth2 ? "AB2-SE-$(filter)" : "RK-SE-$(filter)"
near_global_label(timestepper, filter, fs::ImplicitFreeSurface) = timestepper === :QuasiAdamsBashforth2 ? "AB2-IM" : "RK-IM"

near_global_filename(label, timestepper, filter, fs) = "near_global_" * something(label, near_global_label(timestepper, filter, fs))

function near_global_variants()
    RK3 = RungeKutta3Scheme()
    FB  = ForwardBackwardScheme()
    return [(; label = "AB2-SE",       timestepper = :QuasiAdamsBashforth2, filter = :SM05,   barotropic = FB,  implicit = false),
            (; label = "RK-SE-SM05",   timestepper = :SplitRungeKutta3,     filter = :SM05,   barotropic = RK3, implicit = false),
            (; label = "RK-SE-mu2",    timestepper = :SplitRungeKutta3,     filter = :mu2,    barotropic = RK3, implicit = false),
            (; label = "RK-SE-trig74", timestepper = :SplitRungeKutta3,     filter = :trig74, barotropic = RK3, implicit = false),
            (; label = "RK-SE-trig2",  timestepper = :SplitRungeKutta3,     filter = :trig2,  barotropic = RK3, implicit = false),
            (; label = "RK-IM",        timestepper = :SplitRungeKutta3,     filter = :SM05,   barotropic = RK3, implicit = true)]
end

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
