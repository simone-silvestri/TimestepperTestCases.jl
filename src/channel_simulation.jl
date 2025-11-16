using Oceananigans.Models.VarianceDissipationComputations: VarianceDissipation
using Oceananigans.Operators
using Oceananigans.TimeSteppers: SplitRungeKuttaTimeStepper

const Lx = 1000kilometers # zonal domain length [m]
const Ly = 2000kilometers # meridional domain length [m]

@inline initial_buoyancy(z, p) = p.ΔB * (exp(z / p.h) - exp(-p.Lz / p.h)) / (1 - exp(-p.Lz / p.h))

@inline mask(y, p) = max(0, (y - p.Ly + p.Lsponge) / p.Lsponge)

@inline function buoyancy_relaxation(i, j, k, grid, clock, model_fields, p)
    timescale = p.λt
    z = znode(i, j, k, grid, Center(), Center(), Center())
    y = ynode(i, j, k, grid, Center(), Center(), Center())

    b★ = initial_buoyancy(z, p)
    bᵢ = @inbounds model_fields.b[i, grid.Ny, k]

    return mask(y, p) / timescale * (b★ - bᵢ) 
end

@inline function buoyancy_flux(i, j, grid, clock, model_fields, p)
    y = ynode(i, j, 1, grid, Center(), Center(), Center())
    Q = ifelse(y > p.y_shutoff, zero(grid), p.Qᵇ * cos(3π * y / p.Ly))
    return Q
end

@inline function u_stress(i, j, grid, clock, model_fields, p)
    y = ynode(j, grid, Center())
    return - p.τ * sin(π * y / p.Ly)
end

@inline ϕ²(i, j, k, ϕ) = @inbounds ϕ[i, j, k]^2

default_bottom_height = (x, y) -> y < 1000kilometers ?  5.600000000000001e-15 * y^3 - 8.4e-9 * y^2 - 200 : -3000.0

function default_grid(arch, zstar, bottom_height)

    # number of grid points
    Nx = 200
    Ny = 400
    Nz = 90

    # Ninty levels spacing
    Δz = [10.0 * ones(6)...,
          11.25, 12.625, 14.125, 15.8125, 17.75, 19.9375, 22.375, 25.125, 28.125, 31.625, 35.5, 39.75,
          42.0 * ones(56)...,
          39.75, 35.5, 31.625, 28.125, 25.125, 22.375, 19.9375, 17.75, 15.8125, 14.125, 12.625, 11.25,
          10.0 * ones(4)...]

    z_faces = zeros(Nz+1)
    for k in Nz : -1 : 1
        z_faces[k] = z_faces[k+1] - Δz[Nz - k + 1]
    end

    z_faces = zstar ? MutableVerticalDiscretization(z_faces) : z_faces

    grid = RectilinearGrid(arch,
                        topology = (Periodic, Bounded, Bounded),
                        size = (Nx, Ny, Nz),
                        halo = (6, 6, 6),
                        x = (0, Lx),
                        y = (0, Ly),
                        z = z_faces)

    @info "Built a grid: $grid."

    return isnothing(bottom_height) ? grid : ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))
end

hasclosure(closure, ClosureType) = closure isa ClosureType
hasclosure(closure_tuple::Tuple, ClosureType) = any(hasclosure(c, ClosureType) for c in closure_tuple)

simulation_Δt(::Val{:QuasiAdamsBashforth2}) = 5minutes
simulation_Δt(::Val{:SplitRungeKutta3}) = 10minutes
simulation_Δt(::Val{:SplitRungeKutta6}) = 20minutes

function default_closure()
    mixing_length = CATKEMixingLength()
    turbulent_kinetic_energy_equation = CATKEEquation(Cᵂϵ=1.0)
    return CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation, minimum_tke=1e-7)
end

function run_channel_simulation(; momentum_advection = WENOVectorInvariant(), 
                                    tracer_advection = TimestepperTestCases.tracer_advection, 
                                             closure = default_closure(),
                                               zstar = true,
                                        restart_file = nothing,
                                                arch = CPU(),
                                       bottom_height = nothing,
                                         timestepper = :SplitRungeKutta,
                                                grid = default_grid(arch, zstar, bottom_height),
                                        initial_file = "tIni_80y_90L.bin",
                                            testcase = "0")

    #####
    ##### Boundary conditions
    #####

    α  = 2e-4     # [K⁻¹] thermal expansion coefficient 
    g  = 9.8061   # [m/s²] gravitational constant
    cᵖ = 3994.0   # [J/K]  heat capacity
    ρ  = 999.8    # [kg/m³] reference density

    parameters = (
        Ly  = grid.Ly,
        Lz  = grid.Lz,
        Δy  = grid.Δyᵃᶜᵃ,
        Qᵇ  = 10 / (ρ * cᵖ) * α * g, # buoyancy flux magnitude [m² s⁻³]    
        y_shutoff = 5 / 6 * grid.Ly, # shutoff location for buoyancy flux [m] 
        τ  = 0.1 / ρ,                # surface kinematic wind stress [m² s⁻²]
        μ  = 1.1e-3,                 # bottom drag damping time-scale [-]
        Lsponge = 9 / 10 * Ly,       # sponge region for buoyancy restoring [m]
        ν  = 3e-4,                   # viscosity for "no-slip" lateral boundary conditions
        ΔB = 8 * α * g,              # surface vertical buoyancy gradient [s⁻²]
        H  =  grid.Lz,               # domain depth [m]
        h  = 1000.0,                 # exponential decay scale of stable stratification [m]
        λt = 7.0days                 # relaxation time scale [s]
    )

    # Buoyancy restoring
    buoyancy_restoring = Forcing(buoyancy_relaxation; discrete_form = true, parameters)

    # Top boundary conditions
    u_stress_bc      = FluxBoundaryCondition(u_stress; discrete_form = true, parameters)
    buoyancy_flux_bc = FluxBoundaryCondition(buoyancy_flux, discrete_form = true, parameters = parameters)

    # Drag is added as a forcing to allow both bottom drag _and_ a no-slip BC
    u_bottom_bc = FluxBoundaryCondition(u_quadratic_bottom_drag; discrete_form = true, parameters)
    v_bottom_bc = FluxBoundaryCondition(u_quadratic_bottom_drag; discrete_form = true, parameters)

    u_ibcs = ImmersedBoundaryCondition(bottom = FluxBoundaryCondition(u_immersed_bottom_drag; discrete_form = true, parameters))
    v_ibcs = ImmersedBoundaryCondition(bottom = FluxBoundaryCondition(v_immersed_bottom_drag; discrete_form = true, parameters))
    
    b_bcs = FieldBoundaryConditions(top = buoyancy_flux_bc)
    u_bcs = FieldBoundaryConditions(bottom = u_bottom_bc, immersed = u_ibcs, top = u_stress_bc)
    v_bcs = FieldBoundaryConditions(bottom = v_bottom_bc, immersed = v_ibcs)

    #####
    ##### Coriolis
    #####
    
    actual_Δt = simulation_Δt(Val(timestepper))

    coriolis = BetaPlane(f₀ = -1e-4, β = 1e-11)
    free_surface = SplitExplicitFreeSurface(grid; cfl=0.7, fixed_Δt=actual_Δt)
    tracers = hasclosure(closure, CATKEVerticalDiffusivity) ? (:b, :e) : (:b, )
    
    if closure isa Tuple
       closure = (closure..., VerticalScalarDiffusivity(κ=1e-5, ν=1e-4))
    else
       closure = (closure, VerticalScalarDiffusivity(κ=1e-5,ν=1e-4))
    end

    model = HydrostaticFreeSurfaceModel(; grid,
                                          free_surface,
                                          momentum_advection,
                                          tracer_advection,
                                          buoyancy = BuoyancyTracer(),
                                          coriolis,
                                          closure,
                                          tracers,
                                          timestepper,
                                          forcing = (; b = buoyancy_restoring), 
                                          boundary_conditions = (; b = b_bcs, u = u_bcs, v = v_bcs))

    @info "Built $model."

    #####
    ##### Initial conditions
    #####

    Nx, Ny, Nz = size(grid)

    if  restart_file isa String # Initialize from spinned up solution
        set!(model, restart_file)

    else # resting initial condition
      
      binit = if initial_file isa String
          
          # Initial condition from MITgcm
          Tinit = Array{Float64}(undef, Nx*Ny*Nz)
          read!(initial_file, Tinit)
          Tinit = bswap.(Tinit) |> Array{Float64}
          Tinit = reshape(Tinit, Nx, Ny, Nz)
          binit = reverse(Tinit, dims = 3) .* α .* g

      else
          (x, y, z) -> initial_buoyancy(z, parameters)
      end
          

      set!(model, b = binit) 

    end

    #####
    ##### Simulation building
    #####

    Δt₀ = 1minutes

    # 50 years of simulation
    simulation = Simulation(model; Δt = Δt₀, stop_time = 100days)

    # add progress callback
    wall_clock = [time_ns()]

    function print_progress(sim)
        @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u, w): (%6.3e, %6.3e) m/s, next Δt: %s\n",
            100 * (sim.model.clock.time / sim.stop_time),
            sim.model.clock.iteration,
            prettytime(sim.model.clock.time),
            prettytime(1e-9 * (time_ns() - wall_clock[1])),
            maximum(abs, interior(sim.model.velocities.u)),
            maximum(abs, interior(sim.model.velocities.w)),
            prettytime(sim.Δt))

        wall_clock[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(20))

    # Fuck the spin up!
    if !(restart_file isa String) # Spin up!    
        simulation.output_writers[:first_checkpointer] = Checkpointer(model,
                                                                      schedule = TimeInterval(100days),
                                                                      prefix = "restart" * string(testcase),
                                                                      overwrite_existing = true)
        run!(simulation)

        # Remove first checkpoint
        delete!(simulation.output_writers, :first_checkpointer)

        # Reset time step and simulation time
        model.clock.time = 0
        model.clock.iteration = 0
    end

    simulation.stop_time = 14400days
    simulation.Δt = actual_Δt

    simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                            schedule = TimeInterval(1800days),
                                                            prefix = "channel_checkpoint_" * string(testcase),
                                                            overwrite_existing = true)
    #####
    ##### Diagnostics
    #####
    
    ϵ = Oceananigans.Models.VarianceDissipation(:b, grid)
    simulation.callbacks[:compute_variance] = Callback(ϵ, IterationInterval(1))

    @info "added the tracer variance diagnostic"

    f = Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵ)
    b = model.tracers.b

    VFCC = Oceananigans.AbstractOperations.grid_metric_operation((Face,   Center, Center), Oceananigans.Operators.volume, grid)
    VCFC = Oceananigans.AbstractOperations.grid_metric_operation((Center, Face,   Center), Oceananigans.Operators.volume, grid)
    VCCF = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Face),   Oceananigans.Operators.volume, grid)
    VCCC = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Center), Oceananigans.Operators.volume, grid)

    vol = (; VCCC, VFCC, VCFC, VCCF)

    Gbx = ∂x(b)^2 * VFCC
    Gby = ∂y(b)^2 * VCFC
    Gbz = ∂z(b)^2 * VCCF

    g = (; Gbx, Gby, Gbz)

    snapshot_outputs = merge(model.velocities, model.tracers, f, g, (; η = model.free_surface.η), vol)
    average_outputs  = merge(snapshot_outputs, f, g)

    #####
    ##### Build checkpointer and output writer
    #####

    overwrite_existing = true

    simulation.output_writers[:snapshots] = JLD2Writer(model, snapshot_outputs; 
                                                       schedule = ConsecutiveIterations(TimeInterval(360days)),
                                                       filename = "snapshots_" * string(testcase),
                                                       overwrite_existing)

    simulation.output_writers[:averages] = JLD2Writer(model, average_outputs; 
                                                      schedule = AveragedTimeInterval(5 * 360days),
                                                      filename = "averages_" * string(testcase),
                                                      overwrite_existing)

    @info "Running the simulation..."

    run!(simulation)

    return simulation
end
