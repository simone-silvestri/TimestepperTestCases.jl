using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids
using Oceananigans.ImmersedBoundaries: immersed_inactive_node
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index
using Oceananigans.Operators

@inline function wind_stress(i, j, grid, clock, fields, p) 
    force = clock.time > 3days
    τx = p.τ₀ * sin(0.92 * p.f * clock.time)
    return ifelse(force, τx, zero(grid))
end

idealized_coast_timestep(::Val{:QuasiAdamsBashforth2}) = 60seconds
idealized_coast_timestep(::Val{:SplitRungeKutta3})     = 3minutes

function idealized_coast(timestepper::Symbol; arch = CPU(), forced = true)

    Lx = 97kilometers
    Ly = 97kilometers
    Lz = 103meters

    Nx = 194
    Ny = 194
    Nz = 150

    z_faces = reverse(- [(k / Nz)^(1.25) for k in 0:Nz] .* Lz)
    
    grid = RectilinearGrid(arch; 
                           size = (Nx, Ny, Nz),
                           x = (0, Lx),
                           y = (0, Ly),
                           z = MutableVerticalDiscretization(z_faces),
                           halo = (7, 7, 5),
                           topology = (Periodic, Bounded, Bounded))

    bottom_height(x, y) = - 0.001 * y - 5

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true)

    α = 1.7e-4
    β = 7.6e-4
    f = 1.0e-4

    coriolis = FPlane(; f)

    equation_of_state = LinearEquationOfState(thermal_expansion=α, haline_contraction=β)
    buoyancy = SeawaterBuoyancy(; equation_of_state)

    if forced 
        τ₀ = 0.1 / 1027
    else
        τ₀ = 0.0
    end

    u_top = FluxBoundaryCondition(wind_stress; discrete_form=true, parameters=(τ₀=τ₀, f=f))
    u_bcs = FieldBoundaryConditions(top=u_top)

    model = HydrostaticFreeSurfaceModel(; grid,
                                          coriolis,
                                          timestepper,
                                          tracers = (:T, :S),
                                          buoyancy,
                                          boundary_conditions = (; u=u_bcs,),
                                          free_surface = SplitExplicitFreeSurface(grid; substeps=50),
                                          momentum_advection = WENOVectorInvariant(),
                                          tracer_advection = WENO(order=7))

    N² = 1e-4
    M²(y) = if y > 60kilometers
        0.0
    elseif 50kilometers ≤ y ≤ 60kilometers
        1e-6 - 1e-10 * (y - 50kilometers)
    else
        1e-6
    end
    g  = buoyancy.gravitational_acceleration

    Tᵢ(x, y, z) = 25 + N² / (α * g) * z 
    Sᵢ(x, y, z) = 35 - M²(y) / (β * g) * (y - 60kilometers)
    uᵢ(x, y, z) = y > 60kilometers ? 0.0 : - 1 / f * M²(y) * (z - bottom_height(x, y))

    set!(model, S=Sᵢ, T=Tᵢ, u=uᵢ)
    Δt = idealized_coast_timestep(Val(timestepper))
    simulation = Simulation(model; Δt, stop_time=20days)

    add_callback!(simulation, print_progress,  IterationInterval(100))

    filename = "idealized_coast"
    save_fields_interval = 1hours

    simulation.output_writers[:values] = JLD2Writer(model, merge(model.velocities, model.tracers);
                                                    filename = filename * "_$(string(timestepper))",
                                                    schedule = IterationInterval(300),
                                                    overwrite_existing = true)


    @info "Running the simulation..."

    run!(simulation)

    @info "Simulation completed in " * prettytime(simulation.run_wall_time)

    return simulation
end