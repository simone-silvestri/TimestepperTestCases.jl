using Oceananigans.Models.VarianceDissipationComputations: VarianceDissipation
using Oceananigans.Operators
using Oceananigans.TimeSteppers: SplitRungeKuttaTimeStepper

const Lx = 1000kilometers # zonal domain length [m]
const Ly = 2000kilometers # meridional domain length [m]

"""
    initial_buoyancy(z, p)

Compute the initial buoyancy profile as a function of depth.

$(SIGNATURES)

# Arguments
- `z`: Vertical position [m]
- `p`: Parameters named tuple containing `ΔB` (surface buoyancy gradient), `h` (decay scale), and `Lz` (domain depth)

# Returns
- Initial buoyancy [m/s²] following an exponential profile

The profile follows `ΔB * (exp(z/h) - exp(-Lz/h)) / (1 - exp(-Lz/h))`, which provides
a stable stratification that decays exponentially with depth.
"""
@inline initial_buoyancy(z, p) = p.ΔB * (exp(z / p.h) - exp(-p.Lz / p.h)) / (1 - exp(-p.Lz / p.h))

"""
    mask(y, p)

Compute the sponge layer mask function for buoyancy relaxation.

$(SIGNATURES)

# Arguments
- `y`: Meridional position [m]
- `p`: Parameters named tuple containing `Ly` (domain length) and `Lsponge` (sponge layer width)

# Returns
- Mask value between 0 and 1, increasing linearly from 0 at `y = Ly - Lsponge` to 1 at `y = Ly`

This function creates a sponge layer mask used for buoyancy restoring at the northern boundary.
"""
@inline mask(y, p) = max(0, (y - p.Ly + p.Lsponge) / p.Lsponge)

"""
    buoyancy_relaxation(i, j, k, grid, clock, model_fields, p)

Compute the buoyancy relaxation forcing in the sponge layer.

$(SIGNATURES)

# Arguments
- `i, j, k`: Grid indices
- `grid`: Grid object
- `clock`: Simulation clock
- `model_fields`: Model fields
- `p`: Parameters named tuple containing relaxation parameters

# Returns
- Buoyancy relaxation rate [m/s³] that restores buoyancy toward the initial profile in the sponge layer

This forcing relaxes buoyancy toward the initial exponential profile in the northern sponge region,
with strength controlled by `mask(y, p) / λt` where `λt` is the relaxation timescale.
"""
@inline function buoyancy_relaxation(i, j, k, grid, clock, model_fields, p)
    timescale = p.λt
    z = znode(i, j, k, grid, Center(), Center(), Center())
    y = ynode(i, j, k, grid, Center(), Center(), Center())

    b★ = initial_buoyancy(z, p)
    bᵢ = @inbounds model_fields.b[i, grid.Ny, k]

    return mask(y, p) / timescale * (b★ - bᵢ) 
end

"""
    buoyancy_flux(i, j, grid, clock, model_fields, p)

Compute the surface buoyancy flux forcing.

$(SIGNATURES)

# Arguments
- `i, j`: Grid indices
- `grid`: Grid object
- `clock`: Simulation clock
- `model_fields`: Model fields
- `p`: Parameters named tuple containing `Qᵇ` (flux magnitude), `y_shutoff` (shutoff location), and `Ly` (domain length)

# Returns
- Surface buoyancy flux [m²/s³] following a cosine profile, zero beyond `y_shutoff`

The flux profile is `Qᵇ * cos(3π * y / Ly)` for `y < y_shutoff`, providing differential heating/cooling
that drives meridional circulation.
"""
@inline function buoyancy_flux(i, j, grid, clock, model_fields, p)
    y = ynode(i, j, 1, grid, Center(), Center(), Center())
    Q = ifelse(y > p.y_shutoff, zero(grid), p.Qᵇ * cos(3π * y / p.Ly))
    return Q
end

"""
    u_stress(i, j, grid, clock, model_fields, p)

Compute the zonal wind stress at the surface.

$(SIGNATURES)

# Arguments
- `i, j`: Grid indices
- `grid`: Grid object
- `clock`: Simulation clock
- `model_fields`: Model fields
- `p`: Parameters named tuple containing `τ` (stress magnitude) and `Ly` (domain length)

# Returns
- Zonal wind stress [m²/s²] following a sine profile: `-τ * sin(π * y / Ly)`

This implements a sinusoidal wind stress profile that drives zonal flow, with maximum
stress at the domain center and zero at the boundaries.
"""
@inline function u_stress(i, j, grid, clock, model_fields, p)
    y = ynode(j, grid, Center())
    return - p.τ * sin(π * y / p.Ly)
end

@inline ϕ²(i, j, k, ϕ) = @inbounds ϕ[i, j, k]^2

default_bottom_height = (x, y) -> y < 1000kilometers ?  5.600000000000001e-15 * y^3 - 8.4e-9 * y^2 - 200 : -3000.0

"""
    default_grid(arch, zstar, bottom_height)

Construct the default grid for the channel simulation.

$(SIGNATURES)

# Arguments
- `arch`: Architecture to use (`CPU()` or `GPU()`)
- `zstar`: Whether to use z-star vertical coordinates
- `bottom_height`: Optional bottom height function `(x, y) -> z`, or `nothing` for flat bottom

# Returns
- `RectilinearGrid` or `ImmersedBoundaryGrid` with 200×400×90 grid points

The grid spans 1000 km zonally (periodic), 2000 km meridionally (bounded), and uses
90 non-uniform vertical layers with finer spacing near the surface and bottom.
"""
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

"""
    hasclosure(closure, ClosureType)

Check if a closure or tuple of closures contains a specific closure type.

$(SIGNATURES)

# Arguments
- `closure`: A single closure or tuple of closures
- `ClosureType`: The closure type to check for

# Returns
- `true` if `closure` is of type `ClosureType` or contains it in a tuple, `false` otherwise

This helper function is used to determine if CATKE or other specific closures are present
in the closure configuration.
"""
hasclosure(closure, ClosureType) = closure isa ClosureType
hasclosure(closure_tuple::Tuple, ClosureType) = any(hasclosure(c, ClosureType) for c in closure_tuple)

"""
    simulation_Δt(::Val{timestepper})

Return the recommended time step for the channel simulation given a timestepper.

$(SIGNATURES)

# Arguments
- `timestepper`: Symbol indicating the timestepper (`:QuasiAdamsBashforth2`, `:SplitRungeKutta3`, or `:SplitRungeKutta6`)

# Returns
- Recommended time step [s] for the given timestepper

Time steps are chosen to match computational cost between different timesteppers.
"""
simulation_Δt(::Val{:QuasiAdamsBashforth2}) = 5minutes
simulation_Δt(::Val{:SplitRungeKutta3}) = 10minutes
simulation_Δt(::Val{:SplitRungeKutta6}) = 20minutes

"""
    default_closure()

Return the default turbulence closure for the channel simulation.

$(SIGNATURES)

# Returns
- `CATKEVerticalDiffusivity` closure with default parameters

This closure uses the CATKE (Convective Adjustment and TKE) parameterization for vertical
mixing, which is appropriate for equilibrated channel flows with mesoscale forcing.
"""
function default_closure()
    mixing_length = CATKEMixingLength()
    turbulent_kinetic_energy_equation = CATKEEquation(Cᵂϵ=1.0)
    return CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation, minimum_tke=1e-7)
end

"""
    channel_simulation(; momentum_advection, tracer_advection, closure, zstar, restart_file, arch, bottom_height, timestepper, grid, initial_file, testcase)

Set up and run the idealized re-entrant channel flow simulation.

$(SIGNATURES)

# Keyword Arguments
- `momentum_advection`: Momentum advection scheme (default: `WENOVectorInvariant()`)
- `tracer_advection`: Tracer advection scheme (default: 7th-order WENO)
- `closure`: Turbulence closure (default: `CATKEVerticalDiffusivity`)
- `zstar`: Whether to use z-star vertical coordinates (default: `true`)
- `restart_file`: Optional restart file path (default: `nothing`)
- `arch`: Architecture to run on (default: `CPU()`)
- `bottom_height`: Optional bottom height function (default: `nothing` for flat bottom)
- `timestepper`: Timestepper symbol (default: `:SplitRungeKutta`)
- `grid`: Grid to use (default: `default_grid(...)`)
- `initial_file`: Optional initial condition file (default: `"tIni_80y_90L.bin"`)
- `testcase`: Test case identifier string (default: `"0"`)

# Returns
- `Simulation` object after running to completion

This function sets up the idealized re-entrant channel test case described in the paper.
The configuration consists of a 1000 km × 2000 km × 3 km periodic channel forced by
sinusoidal wind stress and variable surface heat flux. Buoyancy is restored at the northern
boundary to maintain an equilibrated state.

The simulation includes a spin-up phase followed by a long integration (40 years) with
time-averaged outputs. This test case demonstrates how numerical mixing interacts with
explicit physical mixing in an equilibrated configuration, as shown in the paper.
"""
function channel_simulation(; momentum_advection = WENOVectorInvariant(), 
                                tracer_advection = TimestepperTestCases.tracer_advection, 
                                         closure = default_closure(),
                                           zstar = true,
                                    restart_file = nothing,
                                            arch = CPU(),
                                   bottom_height = nothing,
                                     timestepper = :SplitRungeKutta3,
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

    model = HydrostaticFreeSurfaceModel(grid;
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
