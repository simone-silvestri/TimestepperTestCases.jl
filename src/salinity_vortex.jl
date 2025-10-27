using Oceananigans

function salinity_vortex(timestepper::Symbol)
    Nx = 75
    Ny = 75
    Nz = 20

    grid = RectilinearGrid(size=(Nx, Ny, Nz),
                           x = (-15kilometers, 15kilometers),
                           y = (-15kilometers, 15kilometers), 
                           z = MutableVerticalDiscretization((-20meters, 0)),
                           halo = (6, 6, 6),
                           topology = (Bounded, Bounded, Bounded))

    
    f(r) = min(34.85, 1.1 * (r / 3000)^8 + 33.75)
    Sᵢ(x, y, z) = f(sqrt(x^2 + y^2))

    model = HydrostaticFreeSurfaceModel(; grid,
                                         timestepper,
                                         tracer_advection,
                                         momentum_advection = WENOVectorInvariant(),
                                         buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(), constant_temperature=10),
                                         tracers = :S,
                                         coriolis = FPlane(f=1e-4))

    
    set!(model.tracers.S, Sᵢ)

    simulation = Simulation(model, Δt = 60seconds, stop_time = 144hours)

    add_callback!(simulation, print_progress, IterationInterval(100))

    return simulation
end
