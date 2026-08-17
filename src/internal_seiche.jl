using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids

"""
    internal_seiche_parameters()

Return a named tuple containing the default parameters for the internal seiche test case.

$(SIGNATURES)

# Returns
- `Nx`: Number of grid points in the horizontal (x) direction
- `Nz`: Number of grid points in the vertical (z) direction
- `H`: Domain depth [m]
- `L`: Domain length [m]
- `N²`: Background stratification [s⁻²]
- `q`: Vertical mode number to excite (1 is the first baroclinic mode)
- `nwaves`: Number of horizontal wavelengths across the domain
- `u₀`: Velocity amplitude of the initial mode [m/s]

A flat-bottomed, non-rotating, linearly stratified box in which a single vertical normal mode is
released and left to ring. The configuration is the physical realization of the modal analysis of the
manuscript appendix, at its reference stratification `ε = N²H/g = 0.06`.

`Nx` and `nwaves` are set together, and not freely. `μ₀ = k c₀ Δt` is a *grid-scale* quantity: reaching
a large `μ₀` by increasing `Δt` raises the Courant number of every wavenumber on the grid, and the 2Δx
modes -- seeded by round-off -- go unstable long before the excited one does. Stability requires
`α_q μ₀ (k_grid / k_mode) ≲ sqrt(3)`, which caps the number of cells per wavelength. Eight cells per
wavelength keeps the sweep stable to `μ₀ ≈ 7` while holding the discrete phase-speed error to about
10%; see `μ₀_eff` in [`load_internal_seiche`](@ref).

The amplitude is small and the dynamics are linear by construction (no momentum advection, centered
tracer advection, no closure, no drag), so the only dissipation present is that of the time
discretization. This is what makes the test able to resolve per-step damping rates of order `10⁻⁴`.
"""
@inline function internal_seiche_parameters()
    Nx     = 32
    Nz     = 32
    H      = 4kilometers
    L      = 500kilometers
    N²     = 1.4715e-4   # ε = N²H/g = 0.06, the reference stratification of the appendix
    q      = 1           # first baroclinic mode
    nwaves = 4           # 8 cells per wavelength; see the note above
    u₀     = 1e-3        # m/s, small enough that the dynamics stay linear

    return (; Nx, Nz, H, L, N², q, nwaves, u₀)
end

"""
    seiche_vertical_mode(q, param)

Return the exact vertical normal mode `q` of the resting stratification.

$(SIGNATURES)

# Arguments
- `q`: Mode number, `0` for the barotropic mode and `q ≥ 1` for the baroclinic modes
- `param`: Parameters named tuple from [`internal_seiche_parameters`](@ref)

# Returns
Named tuple with

- `xq`: Root of `x tan x = ε`, the phase the mode accumulates over the column
- `αq`: Scaled phase speed, `c_q / sqrt(gH)`
- `cq`: Phase speed [m/s]
- `Aq`: Mean-square normalisation constant
- `φ`: The mode shape `φ(z)`
- `dφ`: Its vertical derivative `dφ/dz`

The modes are kept exact: the barotropic mode is not assumed depth-independent and its speed is not set
to `sqrt(gH)`. That distinction is the whole point of the analysis this case is built to test.
"""
function seiche_vertical_mode(q, param)
    H, N² = param.H, param.N²
    g  = Oceananigans.defaults.gravitational_acceleration
    ε  = N² * H / g

    # Bisect x tan x = ε on the branch containing the qth root.
    lo = q == 0 ? 1e-14 : q * π + 1e-14
    hi = q * π + π/2 - 1e-10
    f(x) = x * tan(x) - ε
    flo = f(lo)
    for _ in 1:200
        mid = (lo + hi) / 2
        fmid = f(mid)
        if flo * fmid ≤ 0
            hi = mid
        else
            lo = mid
            flo = fmid
        end
    end
    xq = (lo + hi) / 2

    αq = sqrt(ε) / xq
    cq = αq * sqrt(g * H)
    Aq = 1 / sqrt(1/2 + sin(2xq) / (4xq))     # mean-square normalisation

    φ(z)  =  Aq * cos(xq * (1 + z / H))
    dφ(z) = -Aq * (xq / H) * sin(xq * (1 + z / H))

    return (; xq, αq, cq, Aq, φ, dφ)
end

"""
    internal_seiche_timestep(μ₀, param)

Return the time step that puts the barotropic Courant number at `μ₀`.

$(SIGNATURES)

# Arguments
- `μ₀`: Target barotropic Courant number `k c₀ Δt`, with `k` the horizontal wavenumber of the excited
  mode and `c₀` the exact barotropic phase speed
- `param`: Parameters named tuple from [`internal_seiche_parameters`](@ref)

# Returns
- Time step [s]

`μ₀` is the coordinate of every result in the appendix, so the sweeps here are specified through it
rather than through `Δt`. The resonance the slow-forcing reconstruction removes sits at `μ₀ ≈ 5.7`.
"""
function internal_seiche_timestep(μ₀, param = internal_seiche_parameters())
    k  = 2π * param.nwaves / param.L
    c₀ = seiche_vertical_mode(0, param).cq
    return μ₀ / (k * c₀)
end

"""
    internal_seiche_grid(; mutable=false)

Construct the grid for the internal seiche test case.

$(SIGNATURES)

# Returns
- `RectilinearGrid`, periodic in `x`, `Flat` in `y`, bounded in `z`, with a flat bottom

The bottom is flat and there is no immersed boundary: the case is deliberately free of every ingredient
that would introduce dissipation of its own, so that the decay measured is the time discretization's.
"""
function internal_seiche_grid(; mutable = false)
    param = internal_seiche_parameters()

    Nx, Nz = param.Nx, param.Nz
    H, L   = param.H, param.L

    z = mutable ? MutableVerticalDiscretization((-H, 0)) : (-H, 0)

    return RectilinearGrid(size = (Nx, Nz), halo = (4, 4),
                           x = (0, L), z = z,
                           topology = (Periodic, Flat, Bounded))
end

"""
    internal_seiche(timestepper::Symbol; μ₀, free_surface, free_surface_name, stop_iteration)

Set up and run the internal seiche test case.

$(SIGNATURES)

# Arguments
- `timestepper`: Symbol indicating the timestepper (`:QuasiAdamsBashforth2` or `:SplitRungeKutta3`)

# Keyword Arguments
- `μ₀`: Barotropic Courant number at which to run (default `5.7`, the resonance)
- `grid`: Grid (default: flat-bottomed box from [`internal_seiche_grid`](@ref))
- `free_surface`: Free surface formulation
- `free_surface_name`: Name used in the output filename
- `stop_iteration`: Number of time steps (default `600`)

# Returns
- `Simulation` object after running to completion

A single vertical mode is released at `t = 0` and left to decay. Because the configuration carries no
physical dissipation, the decay rate of the mode energy measures the amplification factor of the coupled
barotropic--baroclinic scheme directly, one factor of `exp(2λ)` per step.

Two things are read off a sweep over `μ₀`. At small `μ₀` the error against the exact modal solution gives
the *order* of the scheme, which distinguishes the forward--backward substep (first order) from the
three-stage one (second). Near `μ₀ ≈ 5.7` the measured damping shows the resonance of the frozen slow
forcing, which loses about a third of its damping there and which the stage-value reconstruction removes.

The run is deliberately tiny (32 × 32) and cheap enough to sweep densely in `μ₀`. The excited wave is
only a few cells across, which is the regime the appendix's `k = 1/Δx` describes.
"""
function internal_seiche(timestepper::Symbol;
                         μ₀ = 5.7,
                         grid = internal_seiche_grid(),
                         free_surface = SplitExplicitFreeSurface(grid; substeps = 48,
                                                                 timestepper = RungeKutta3Scheme(),
                                                                 averaging_kernel = OptimizedAsymmetricAveragingKernel()),
                         free_surface_name = default_free_surface_name(free_surface),
                         stop_iteration = 600,
                         save_interval = 1)

    param = internal_seiche_parameters()
    mode  = seiche_vertical_mode(param.q, param)
    k     = 2π * param.nwaves / param.L

    # Linear by construction: no momentum advection, centered tracer advection, no closure, no rotation.
    model = HydrostaticFreeSurfaceModel(grid;
                                        buoyancy = BuoyancyTracer(),
                                        tracers = :b,
                                        momentum_advection = nothing,
                                        tracer_advection = Centered(),
                                        free_surface,
                                        timestepper)

    # Mode q released at rest phase: p = A φ(z) cos(kx), u = -(A/c) φ(z) cos(kx), b' = A φ'(z) cos(kx).
    A = param.u₀ * mode.cq / mode.Aq

    uᵢ(x, z) = -(A / mode.cq) * mode.φ(z) * cos(k * x)
    bᵢ(x, z) = param.N² * z + A * mode.dφ(z) * cos(k * x)
    # The free surface lives at (Center, Center, Face), so `set!` hands it (x, z).
    ηᵢ(x, z) = (A / Oceananigans.defaults.gravitational_acceleration) * mode.φ(0) * cos(k * x)

    set!(model, u = uᵢ, b = bᵢ)
    set!(model.free_surface.displacement, ηᵢ)

    Δt = internal_seiche_timestep(μ₀, param)
    simulation = Simulation(model; Δt, stop_iteration)

    add_callback!(simulation, print_progress, IterationInterval(200))

    u, v, w = model.velocities
    b = model.tracers.b
    η = model.free_surface.displacement

    # Perturbation energy. The background is removed by subtracting the horizontal mean, which is exact
    # here because the mode is a single cosine and integrates to zero over the domain.
    b̄  = Field(Average(b, dims = 1))
    b′ = b - b̄

    KE   = Field(Average(u^2 / 2))
    APE  = Field(Average(b′^2 / (2 * param.N²)))
    ηvar = Field(Average(η^2))

    filename = "internal_seiche_$(string(timestepper))_$(free_surface_name)_mu$(replace(string(μ₀), "." => "p"))"

    simulation.output_writers[:energy] = JLD2Writer(model, (; KE, APE, ηvar);
                                                    filename = filename * "_energy",
                                                    schedule = IterationInterval(save_interval),
                                                    overwrite_existing = true)

    simulation.output_writers[:fields] = JLD2Writer(model, (; u, w, b, η);
                                                    filename,
                                                    schedule = IterationInterval(max(1, stop_iteration ÷ 20)),
                                                    overwrite_existing = true)

    run!(simulation)

    return simulation
end
