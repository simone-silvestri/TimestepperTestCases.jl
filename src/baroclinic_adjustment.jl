using Oceananigans
using Oceananigans.Units
using Random

baroclinic_timestep(::Val{:QuasiAdamsBashforth2}) = 20minutes
baroclinic_timestep(::Val{:SplitRungeKutta3})     = 60minutes
baroclinic_timestep(::Val{:SplitRungeKutta6})     = 120minutes

function baroclinic_adjustment(timestepper::Symbol)
    Lx = 1000kilometers 
    Ly = 1000kilometers 
    Lz = 1kilometers    

    grid = RectilinearGrid(size = (48, 48, 8),
                           x = (-Lx/2, Lx/2),
                           y = (-Ly/2, Ly/2),
                           z = (-Lz, 0), #MutableVerticalDiscretization((-Lz, 0)),
                           halo = (6, 6, 6),
                           topology = (Periodic, Bounded, Bounded))

    # bottom_height = (- Lz - 10) .* ones(48, 50)
    # bottom_height[:, 1] .= 0
    # bottom_height[:, end] .= 0

    # grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

    model = HydrostaticFreeSurfaceModel(; grid,
                                          coriolis = BetaPlane(latitude = -45),
                                          buoyancy = BuoyancyTracer(),
                                          timestepper,
                                          tracers = :b,
                                          free_surface = SplitExplicitFreeSurface(grid; substeps=50),
                                        #   closure = HorizontalScalarBiharmonicDiffusivity(ν=2e10),
                                          momentum_advection = WENO(),
                                          tracer_advection = WENO(order=7))

    ramp(y, Δy) = min(max(0, y/Δy + 1/2), 1)

    N² = 1e-5
    M² = 1e-7

    Random.seed!(1234)
    Δy = 100kilometers 
    Δb = Δy * M²       
    ϵb = 1e-5

    bᵢ(x, y, z) = N² * z + Δb * ramp(y, Δy) + ϵb * randn()

    Δt = baroclinic_timestep(Val(timestepper))
    set!(model, b=bᵢ)
    
    simulation = Simulation(model, Δt=Δt, stop_time=200days)
    add_callback!(simulation, print_progress, IterationInterval(100))

    u, v, _ = model.velocities
    b = model.tracers.b

    ζ  = ∂x(v) - ∂y(u)
    bx = ∂x(b)^2
    by = ∂y(b)^2
    bz = ∂z(b)^2

    filename = "baroclinic_adjustment"
    save_fields_interval = 0.5day

    ϵ = Oceananigans.Models.VarianceDissipationComputations.VarianceDissipation(:b, grid)
    f = Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵ)
    add_callback!(simulation, ϵ, IterationInterval(1))

    simulation.output_writers[:values] = JLD2Writer(model, merge(model.velocities, model.tracers, (; ζ, bx, by, bz), f);
                                                    filename = filename * "_$(string(timestepper))",
                                                    schedule = TimeInterval(save_fields_interval),
                                                    overwrite_existing = true)

    @info "Running the simulation..."

    run!(simulation)

    @info "Simulation completed in " * prettytime(simulation.run_wall_time)

    return simulation
end
