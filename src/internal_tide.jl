using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids

"""
    internal_tide_parameters()

Return a named tuple containing the default parameters for the internal tide test case.

$(SIGNATURES)

# Returns
- `Nx`: Number of grid points in the horizontal (x) direction
- `Nz`: Number of grid points in the vertical (z) direction
- `H`: Domain depth [m]
- `L`: Domain half-length [m]
- `h₀`: Seamount height [m]
- `width`: Seamount width [m]
- `T₂`: Tidal period [s]
- `ω₂`: Tidal frequency [rad/s]
- `ϵ`: Excursion parameter (dimensionless)
- `U₂`: Characteristic tidal velocity [m/s]
- `f`: Coriolis parameter [s⁻¹]
- `A₂`: Tidal forcing amplitude [m/s²]
- `Nᵢ²`: Initial stratification [s⁻²]

These parameters define a 2-kilometer deep domain with a Gaussian seamount that interacts
with an oscillatory tidal forcing, as described in the paper.
"""
@inline function internal_tide_parameters() 
    Nx    = 256
    Nz    = 128
    H     = 2kilometers
    L     = 1000kilometers
    h₀    = 250meters
    width = 20kilometers
    T₂    = 12.421hours
    ω₂    = 2π / T₂ # radians/sec
    ϵ     = 0.1 # excursion parameter
    U₂    = ϵ * ω₂ * width
    f     = -0.000103126 # coriolis parameter
    A₂    = U₂ * (ω₂^2 - f^2) / ω₂
    Nᵢ²   = 1e-4 # initial stratification (s⁻²)

    return (; Nx, Nz, H, L, h₀, width, T₂, ω₂, ϵ, U₂, f, A₂, Nᵢ²)
end

"""
    tidal_forcing(x, z, t, p)

Compute the tidal forcing amplitude at position `(x, z)` and time `t`.

$(SIGNATURES)

# Arguments
- `x`: Horizontal position [m]
- `z`: Vertical position [m]
- `t`: Time [s]
- `p`: Parameters named tuple containing `A₂` (forcing amplitude) and `ω₂` (tidal frequency)

# Returns
- Forcing amplitude [m/s²] as a function of time only: `A₂ * sin(ω₂ * t)`

This function implements a simple oscillatory tidal forcing that varies sinusoidally in time
with amplitude `A₂` and frequency `ω₂`.
"""
@inline tidal_forcing(x, z, t, p) = p.A₂ * sin(p.ω₂ * t)

"""
    internal_tide_timestep(::Val{timestepper})

Return the recommended time step for the internal tide test case given a timestepper.

$(SIGNATURES)

# Arguments
- `timestepper`: Symbol indicating the timestepper (`:QuasiAdamsBashforth2`, `:SplitRungeKutta3`, or `:SplitRungeKutta4`)

# Returns
- Recommended time step [s] for the given timestepper

The time steps are chosen to match computational cost between AB2 and RK schemes while
maintaining stability, as described in the paper.
"""
internal_tide_timestep(::Val{:QuasiAdamsBashforth2}) =  5minutes
internal_tide_timestep(::Val{:SplitRungeKutta3})     = 10minutes
internal_tide_timestep(::Val{:SplitRungeKutta4})     = 10minutes

@kernel function _compute_dissipation!(Δtc², c⁻, c, Δt)
    i, j, k = @index(Global, NTuple)
    @inbounds begin
        Δtc²[i, j, k] = (c[i, j, k]^2 - c⁻[i, j, k]^2) / Δt
        c⁻[i, j, k]   = c[i, j, k]
    end
end

"""
    compute_tracer_dissipation!(sim)

Compute the time rate of change of tracer variance (dissipation) for tracers `b` and `c`.

$(SIGNATURES)

# Arguments
- `sim`: The `Simulation` object containing the model with tracers `b` and `c`

# Returns
- `nothing` (modifies `sim.model.auxiliary_fields` in place)

This function computes the dissipation rate of tracer variance, defined as
`Δtc² = (c² - c⁻²) / Δt` where `c⁻` is the tracer value from the previous time step.
The results are stored in `sim.model.auxiliary_fields.Δtc²` and `sim.model.auxiliary_fields.Δtb²`.

This diagnostic is used to quantify numerical mixing introduced by the time discretization,
as described in the paper's appendix.
"""
function compute_tracer_dissipation!(sim)
    c    = sim.model.tracers.c
    c⁻   = sim.model.auxiliary_fields.c⁻
    Δtc² = sim.model.auxiliary_fields.Δtc²
    Oceananigans.Utils.launch!(CPU(), sim.model.grid, :xyz,
                               _compute_dissipation!,
                               Δtc², c⁻, c, sim.Δt)

    b    = sim.model.tracers.b
    b⁻   = sim.model.auxiliary_fields.b⁻
    Δtb² = sim.model.auxiliary_fields.Δtb²
    Oceananigans.Utils.launch!(CPU(), sim.model.grid, :xyz,
                               _compute_dissipation!,
                               Δtb², b⁻, b, sim.Δt)

    return nothing
end

"""
    internal_tide_grid()

Construct the grid for the internal tide test case with a Gaussian seamount.

$(SIGNATURES)

# Returns
- `ImmersedBoundaryGrid` with periodic horizontal boundaries and a Gaussian seamount bottom

The grid spans from `-L` to `L` horizontally (periodic) and from `-H` to `0` vertically (bounded).
A Gaussian seamount of height `h₀` and width `width` is imposed using an immersed boundary method.
The grid parameters are taken from `internal_tide_parameters()`.
"""
function internal_tide_grid()
    param = internal_tide_parameters()

    Nx, Nz    = param.Nx, param.Nz
    h₀, width = param.h₀, param.width
    H, L      = param.H, param.L

    underlying_grid = RectilinearGrid(size = (Nx, Nz), halo = (6, 6),
                                    x = (-L, L), z = (-H, 0), # MutableVerticalDiscretization((-H, 0)),
                                    topology = (Periodic, Flat, Bounded))

    hill(x)   =   h₀ * exp(-x^2 / 2width^2)
    bottom(x) = - H + hill(x)

    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom); active_cells_map = true)
   
    return grid
end

"""
    internal_tide(timestepper::Symbol; free_surface, tracer_advection)

Set up and run the internal tide test case simulation.

$(SIGNATURES)

# Arguments
- `timestepper`: Symbol indicating the timestepper (`:QuasiAdamsBashforth2`, `:SplitRungeKutta3`, or `:SplitRungeKutta4`)

# Keyword Arguments
- `free_surface`: Free surface formulation (default: `SplitExplicitFreeSurface` with 60 substeps)
- `tracer_advection`: Tracer advection scheme (default: 7th-order WENO)

# Returns
- `Simulation` object after running to completion

This function sets up the internal tide test case described in the paper, which simulates
tidal flow over a Gaussian seamount. The domain is initially stratified and forced by an
oscillatory tidal forcing. The simulation runs for 40 days and outputs velocity, buoyancy,
tracer fields, and dissipation diagnostics.

The test case isolates the role of time discretization in numerical mixing, as spatial
advection plays a secondary role in this mostly linear configuration.
"""
function internal_tide(timestepper::Symbol; 
                       free_surface=SplitExplicitFreeSurface(internal_tide_grid(); substeps=60),
                       tracer_advection=TimestepperTestCases.tracer_advection)
    
    grid  = internal_tide_grid()
    param = internal_tide_parameters()

    coriolis  = FPlane(f = param.f)
    u_forcing = Forcing(tidal_forcing, parameters=param)

    c⁻    = CenterField(grid)
    Δtc²  = CenterField(grid)
    b⁻    = CenterField(grid)
    Δtb²  = CenterField(grid)

    model = HydrostaticFreeSurfaceModel(; grid, coriolis = nothing,
                                          buoyancy = BuoyancyTracer(),
                                          tracers = (:b, :c),
                                          momentum_advection = WENO(),
                                          tracer_advection = tracer_advection,
                                          free_surface,
                                          timestepper,
                                          forcing = (; u = u_forcing),
                                          auxiliary_fields=(; Δtc², c⁻, Δtb², b⁻))

    bᵢ(x, z) = param.Nᵢ² * z
    cᵢ(x, z) = exp( - (z + 1kilometers)^2 / (2 * (25meters)^2))
    set!(model, u=param.U₂, b=bᵢ)

    Δt = internal_tide_timestep(Val(timestepper)) 
    stop_time = 40days
    simulation = Simulation(model; Δt, stop_time)

    ϵb = Oceananigans.Models.VarianceDissipationComputations.VarianceDissipation(:b, grid)
    ϵc = Oceananigans.Models.VarianceDissipationComputations.VarianceDissipation(:c, grid)

    # Adding the variance dissipation
    add_callback!(simulation, ϵb, IterationInterval(1))
    add_callback!(simulation, ϵc, IterationInterval(1))
    add_callback!(simulation, compute_tracer_dissipation!, IterationInterval(1))

    wall_clock = Ref(time_ns())

    add_callback!(simulation, print_progress, IterationInterval(200))

    u, v, w = model.velocities
    b  = model.tracers.b
    c  = model.tracers.c
    η  = model.free_surface.η
    U  = Field(Average(u))
    u′ = u - U
    N² = ∂z(b)

    Gbx = ∂x(b)^2
    Gbz = ∂z(b)^2
    Gcx = ∂x(c)^2
    Gcz = ∂z(c)^2

    g = (; Gbx, Gbz, Gcx, Gcz)

    if free_surface isa SplitExplicitFreeSurface
        fsname = "split_free_surface_DFB"
    else
        fsname = "implicit_free_surface"
    end

    if tracer_advection != TimestepperTestCases.tracer_advection
        fsname *= "_$(typeof(tracer_advection).name.name)"
    end

    filename = "internal_tide_$(string(timestepper))_$(fsname)"
    save_fields_interval = 1hours
    
    f = merge(Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵb),
              Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵc))

    VFC = Oceananigans.AbstractOperations.grid_metric_operation((Face,   Center, Center), Oceananigans.Operators.volume, grid)
    VCF = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Face),   Oceananigans.Operators.volume, grid)
    VCC = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Center), Oceananigans.Operators.volume, grid)

    V  = (; VFC, VCF, VCC)
    Δ² = (; Δtc² = model.auxiliary_fields.Δtc², Δtb² = model.auxiliary_fields.Δtb²)

    outputs = ( u = u * VFC,
                w = w * VCF,
                b = b * VCC,
                c = c * VCC,
                η = η,
                Δtc² = Δ².Δtc² * VCC,
                Δtb² = Δ².Δtb² * VCC,
                Gbx = Gbx * VFC,
                Gbz = Gbz * VCF,
                Gcx = Gcx * VFC,
                Gcz = Gcz * VCF,
                Abx = f.Abx,
                Abz = f.Abz,
                Acx = f.Acx,
                Acz = f.Acz)

    simulation.output_writers[:fields] = JLD2Writer(model, outputs; 
                                                    filename,
                                                    schedule = TimeInterval(save_fields_interval),
                                                    overwrite_existing = true)

    run!(simulation)

    return simulation
end
