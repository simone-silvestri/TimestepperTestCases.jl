using Oceananigans
using Oceananigans.Utils: ConsecutiveIterations

plume_spreading_timestep(::Val{:QuasiAdamsBashforth2}) = 60seconds
plume_spreading_timestep(::Val{:SplitRungeKutta3})     = 3minutes
plume_spreading_timestep(::Val{:SplitRungeKutta6})     = 9minutes

function plume_spreading_grid(arch)

    Ly = 100kilometers
    Lz = 30meters
    Ny = 200

    # Non-uniform grid in x
    x_faces = [0.0]
    for i in 1:70
        push!(x_faces, x_faces[end] + 1kilometer)
    end

    a = -9.978070175438596

    for i in 71:71+94
        dx = a * (i - 70) + 1kilometer
        push!(x_faces, x_faces[end] + dx)
    end

    for i in 71+95:71+95+19
        push!(x_faces, x_faces[end] + 50meters)
    end

    for i in 71+95+20:71+95+20+94
        dx = 50meters + a * (70 + 95 + 20 - i) 
        push!(x_faces, x_faces[end] + dx)
    end

    for i in 1:110
        push!(x_faces, x_faces[end] + 1kilometer)
    end

    Nx = length(x_faces) - 1

    Nz = 40
    z_faces = - reverse([(k / Nz)^2 for k in 0:Nz]) .* Lz

    grid = RectilinearGrid(arch; 
                           size = (Nx, Ny, Nz),
                           x = x_faces,
                           y = (0, Ly),
                           z = MutableVerticalDiscretization(z_faces),
                           halo = (7, 7, 5),
                           topology = (Bounded, Bounded, Bounded))

    bottom_height(x, y) = if y < 10kilometers
        if 119.5kilometers ≤ x ≤ 120.5kilometers 
            -10meters
        else
            0meters
        end
    elseif y < 20kilometers
        - 10meters - 20meters * (y - 10kilometers) / 10kilometers
    else
        -30meters
    end

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true)
    
    return grid
end

@inline u_restoring(i, j, k, grid, clock, fields, p) = @inbounds ifelse(j ≤ 5, 1 / p.λ * (p.u₀ - fields.u[i, j, k]), zero(u)) * min(1, (clock.time - 1hours) / 1hours)
@inline S_restoring(i, j, k, grid, clock, fields, p) = @inbounds ifelse(j ≤ 5, 1 / p.λ * (0    - fields.S[i, j, k]), zero(S)) * min(1, (clock.time - 1hours) / 1hours)

@inline u_open(i, k, grid, clock, fields, p) = p.u₀ * min(1, (clock.time - 1hours) / 1hours)
@inline S_open(i, k, grid, clock, fields, p) = p.S₀ * max(0, (1hours - clock.time) / 1hours)

function plume_spreading_model(timestepper::Symbol; arch = CPU())

    grid = plume_spreading_grid(arch)
    coriolis = FPlane(f = 1.2e-4)

    parameters = (; λ = 1 / 20minutes, u₀ = 0.3, S₀ = 30.0)
    u_rest = Forcing(u_restoring; discrete_form=true, parameters) 
    S_rest = Forcing(S_restoring; discrete_form=true, parameters)

    u_in =  OpenBoundaryCondition(0.3) # u_open; discrete_form=true, parameters)
    S_in = ValueBoundaryCondition(0.0) # S_open; discrete_form=true, parameters)
    
    u_bcs = FieldBoundaryConditions(south=u_in)
    S_bcs = FieldBoundaryConditions(south=S_in)

    equation_of_state = LinearEquationOfState(thermal_expansion = 0.0)
    buoyancy = SeawaterBuoyancy(; equation_of_state, constant_temperature = 20.0)

    model = HydrostaticFreeSurfaceModel(; grid,
                                          coriolis,
                                          timestepper,
                                          tracers = :S,
                                          forcing = (; u = u_rest, S = S_rest),
                                          buoyancy,
                                          boundary_conditions = (; u = u_bcs, S = S_bcs),
                                          free_surface = SplitExplicitFreeSurface(grid; substeps=10),
                                          momentum_advection = WENOVectorInvariant(),
                                          tracer_advection = WENO(order=7))

    set!(model, S=30.0)
    
    Δt = plume_spreading_timestep(Val(timestepper))

    simulation = Simulation(model; Δt, stop_time=35hours)
    add_callback!(simulation, print_progress, IterationInterval(10))

    filename = "plume_spreading"
    save_fields_interval = 10minutes

    ϵ = Oceananigans.Models.VarianceDissipationComputations.VarianceDissipation(:S, grid)
    f = Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵ)
    add_callback!(simulation, ϵ, ConsecutiveIterations(IterationInterval(300)))

    simulation.output_writers[:values] = JLD2Writer(model, merge(model.velocities, model.tracers, f);
                                                    filename = filename * "_$(string(timestepper))",
                                                    schedule = IterationInterval(300),
                                                    overwrite_existing = true)

    @info "Running the simulation..."

    # run!(simulation)

    @info "Simulation completed in " * prettytime(simulation.run_wall_time)

    return simulation
end