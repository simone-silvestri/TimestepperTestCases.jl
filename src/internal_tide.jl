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
internal_tide_timestep(::Val{scheme}) where scheme =
    baroclinic_timestep(scheme; internal_tide_stability_parameters()...)

"""
    internal_tide_stability_parameters()

Stratification, depth and grid spacing that set the first baroclinic Courant number of this case, and with it
the time step of every scheme through [`baroclinic_timestep`](@ref).
"""
function internal_tide_stability_parameters()
    p = internal_tide_parameters()
    return (N² = p.Nᵢ², H = p.H, Δx = 2p.L / p.Nx, safety = 1)
end

"""
    internal_tide_substeps(barotropic_scheme, scheme; averaging_kernel)

Barotropic substep count for this case, fixed so that the substep Courant number `c₀ k Δτ` sits at 70% of the
limit of `barotropic_scheme` -- `√3` for the three-stage Runge-Kutta substep, `1.8` for forward-backward. The
count therefore depends on the substep integrator as well as on the baroclinic time step.
"""
function internal_tide_substeps(barotropic_scheme, scheme = :SplitRungeKutta3;
                                averaging_kernel = OptimizedAsymmetricAveragingKernel())
    p  = internal_tide_stability_parameters()
    Δt = internal_tide_timestep(Val(scheme))
    return barotropic_substeps(barotropic_scheme; p.H, p.Δx, Δt, averaging_kernel)
end

@kernel function _compute_dissipation!(Δtσc², σc²⁻, c, grid, Δt)
    i, j, k = @index(Global, NTuple)
    @inbounds begin
        σc² = volume(i, j, k, grid, Center(), Center(), Center()) * c[i, j, k]^2
        Δtσc²[i, j, k] = (σc² - σc²⁻[i, j, k]) / Δt
        σc²⁻[i, j, k]  = σc²
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
    grid = sim.model.grid

    c    = sim.model.tracers.c
    c⁻   = sim.model.auxiliary_fields.c⁻      # holds the previous thickness-weighted variance σc²⁻
    Δtc² = sim.model.auxiliary_fields.Δtc²
    Oceananigans.Utils.launch!(CPU(), grid, :xyz,
                               _compute_dissipation!,
                               Δtc², c⁻, c, grid, sim.Δt)

    b    = sim.model.tracers.b
    b⁻   = sim.model.auxiliary_fields.b⁻      # holds the previous thickness-weighted variance σb²⁻
    Δtb² = sim.model.auxiliary_fields.Δtb²
    Oceananigans.Utils.launch!(CPU(), grid, :xyz,
                               _compute_dissipation!,
                               Δtb², b⁻, b, grid, sim.Δt)

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
function internal_tide_grid(; mutable=true)
    param = internal_tide_parameters()

    Nx, Nz    = param.Nx, param.Nz
    h₀, width = param.h₀, param.width
    H, L      = param.H, param.L

    z = mutable ? MutableVerticalDiscretization((-H, 0)) : (-H, 0)

    underlying_grid = RectilinearGrid(size = (Nx, Nz), halo = (6, 6),
                                    x = (-L, L), z = z,
                                    topology = (Periodic, Flat, Bounded))

    hill(x)   =   h₀ * exp(-x^2 / 2width^2)
    bottom(x) = - H + hill(x)

    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom))
   
    return grid
end

default_free_surface_name(::SplitExplicitFreeSurface) = "split_free_surface"
default_free_surface_name(::ImplicitFreeSurface) = "implicit_free_surface"

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
                       grid = internal_tide_grid(),
                       barotropic_timestepper=ForwardBackwardScheme(),
                       averaging_kernel=OptimizedAsymmetricAveragingKernel(),
                       free_surface=SplitExplicitFreeSurface(grid; averaging_kernel,
                                                             timestepper=barotropic_timestepper,
                                                             substeps=internal_tide_substeps(barotropic_timestepper,
                                                                                             timestepper;
                                                                                             averaging_kernel)),
                       free_surface_name=default_free_surface_name(free_surface),
                       tracer_advection=TimestepperTestCases.tracer_advection)

    param = internal_tide_parameters()

    coriolis  = FPlane(f = param.f)
    u_forcing = Forcing(tidal_forcing, parameters=param)

    c⁻    = CenterField(grid)
    Δtc²  = CenterField(grid)
    b⁻    = CenterField(grid)
    Δtb²  = CenterField(grid)

    model = HydrostaticFreeSurfaceModel(grid; 
                                        coriolis,
                                        buoyancy = BuoyancyTracer(),
                                        tracers = (:b, :c),
                                        momentum_advection = WENO(minimum_buffer_upwind_order=1),
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
    η  = model.free_surface.displacement
    U  = Field(Average(u))
    u′ = u - U
    N² = ∂z(b)

    Gbx = ∂x(b)^2
    Gbz = ∂z(b)^2
    Gcx = ∂x(c)^2
    Gcz = ∂z(c)^2

    g = (; Gbx, Gbz, Gcx, Gcz)

    # The advection scheme is appended only when the caller accepted the default name, which is what keeps a
    # non-default-advection run from overwriting the default one. An explicitly supplied name already
    # distinguishes the run, so appending to it would only make the filename disagree with the case label.
    fsname = free_surface_name
    if free_surface_name == default_free_surface_name(free_surface) &&
       tracer_advection != TimestepperTestCases.tracer_advection
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
                Δtc² = Δ².Δtc²,
                Δtb² = Δ².Δtb²,
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

"""
    internal_tide(d::Discretization; grid = internal_tide_grid())

Run the internal tide case with the discretization `d` of Table 1. The time step and the barotropic substep
count are derived from `d` and from `internal_tide_stability_parameters()`.
"""
function internal_tide(d::Discretization; grid = internal_tide_grid())
    _, free_surface = timestep_and_free_surface(d, grid, internal_tide_stability_parameters())
    return internal_tide(d.timestepper; grid, free_surface, free_surface_name = d.label, tracer_advection = d.tracer_advection)
end
