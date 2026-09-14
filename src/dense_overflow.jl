using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids
using Oceananigans.BuoyancyFormulations: LinearEquationOfState

"""
    dense_overflow_parameters()

Geometry, water masses and viscosity of the overflow benchmark of Ilıcak, Adcroft, Griffies and Hallberg
(2012), as the MPAS-Ocean registry states it.

$(SIGNATURES)

# Returns
- `Lx`, `H`: length and depth of the domain [m]
- `Hˢ`, `xˢ`, `Lˢ`: depth of the shelf, position of the shelf break and width of the slope [m]
- `Lᵖ`: width of the dense plug on the shelf [m]
- `Nx`, `Nz`: 1 km by 20 m
- `Tᵖ`, `Tᵃ`, `S`: plug and ambient temperature [°C] and the uniform salinity [psu]
- `α`, `β`: expansion coefficients of the linear equation of state [°C⁻¹, psu⁻¹]
- `U`, `W`: horizontal and vertical speed the time step is sized on, above what the plume reaches [m s⁻¹]
- `ν`: lateral Laplacian viscosity [m² s⁻¹], the low end of the benchmark sweep, the WENO momentum
  advection supplying the grid-scale dissipation the reference value stands in for
- `Cᴰ`: quadratic bottom drag coefficient
"""
@inline function dense_overflow_parameters()
    Lx = 200kilometers
    H  = 2000meters
    Hˢ = 500meters
    xˢ = 40kilometers
    Lˢ = 7kilometers
    Lᵖ = 20kilometers
    Nx = 100
    Nz = 50
    Tᵖ = 10
    Tᵃ = 20
    S  = 35
    α  = 2e-4
    β  = 8e-4
    U  = 1.5
    W  = 0.5
    ν  = 100.0
    Cᴰ = 1e-3

    return (; Lx, H, Hˢ, xˢ, Lˢ, Lᵖ, Nx, Nz, Tᵖ, Tᵃ, S, α, β, U, W, ν, Cᴰ)
end

"""
    dense_overflow_bottom_height(x, p = dense_overflow_parameters())

Sea floor [m], a shelf at `-Hˢ` joined to the basin at `-H` by a tanh of width `Lˢ` centered on `xˢ`.
"""
@inline dense_overflow_bottom_height(x, p = dense_overflow_parameters()) = -(p.H + (p.Hˢ - p.H) / 2 * (1 + tanh((p.xˢ - x) / p.Lˢ)))

"""
    dense_overflow_stability_parameters(p = dense_overflow_parameters())

Stratification, depth, spacing and speed that set the time step of this case.

$(SIGNATURES)

`N²` is the equivalent uniform stratification of the two-layer front, `(π c / H)²` for the interfacial speed
`c = √(g′ Hˢ (H - Hˢ) / H)`, so that `first_baroclinic_speed` returns `c`. The speed is `U + W Δx/Δz`: on a
1 km by 20 m grid the vertical Courant number binds, so the vertical velocity enters scaled to the horizontal
spacing the limit is written on.
"""
function dense_overflow_stability_parameters(p = dense_overflow_parameters())
    Δx = p.Lx / p.Nx
    Δz = p.H / p.Nz
    g′ = Oceananigans.defaults.gravitational_acceleration * p.α * (p.Tᵃ - p.Tᵖ)
    c  = sqrt(g′ * p.Hˢ * (p.H - p.Hˢ) / p.H)
    U  = p.U + p.W * Δx / Δz

    return (; N² = (π * c / p.H)^2, H = p.H, Δx, U, horizontal_dimensions = 1)
end

"""
    dense_overflow_timestep(::Val{timestepper})

Time step [s] of this case, each scheme at the same fraction of its own stability limit.
"""
dense_overflow_timestep(::Val{scheme}) where scheme = baroclinic_timestep(scheme; dense_overflow_stability_parameters()...)

"""
    dense_overflow_substeps(barotropic_scheme, scheme; averaging_kernel)

Barotropic substep count that puts the substep Courant number at the `barotropic_cfl` of `barotropic_scheme`.
The domain is `Flat` in `y`, so the two-point wavenumber of the sub-cycle is `2/Δx`.
"""
function dense_overflow_substeps(barotropic_scheme, scheme = :SplitRungeKutta3;
                                 averaging_kernel = OptimizedAsymmetricAveragingKernel())
    p = dense_overflow_stability_parameters()
    return barotropic_substeps(barotropic_scheme; p.H, p.Δx, averaging_kernel, wavenumber = 2 / p.Δx,
                               Δt = dense_overflow_timestep(Val(scheme)))
end

"""
    dense_overflow_grid(arch = CPU(); p = dense_overflow_parameters())

The `x`-`z` grid of this case, `Bounded` in both directions and `Flat` in `y`, over the shelf and slope of
`dense_overflow_bottom_height`.
"""
function dense_overflow_grid(arch = CPU(); p = dense_overflow_parameters())
    grid = RectilinearGrid(arch; size = (p.Nx, p.Nz), halo = (7, 5),
                           x = (0, p.Lx), z = MutableVerticalDiscretization((-p.H, 0)),
                           topology = (Bounded, Flat, Bounded))

    bottom_height(x) = dense_overflow_bottom_height(x, p)

    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true)
end

"""
    dense_overflow(timestepper::Symbol; kw...)

Set up and run the dense overflow test case.

$(SIGNATURES)

# Arguments
- `timestepper`: the baroclinic composition, as `HydrostaticFreeSurfaceModel` takes it

# Keyword arguments
- `arch`: architecture to run on
- `grid`: the grid, from `dense_overflow_grid`
- `free_surface`, `averaging_kernel`, `barotropic_timestepper`, `slow_forcing`: the free-surface treatment
- `boundary_scheme`, `tracer_boundary_scheme`, `momentum_boundary_scheme`: the reconstruction the WENO buffer
  chain terminates in, as in `idealized_coast`
- `tracer_advection`, `momentum_advection`: built from the boundary schemes where not given
- `Δt`: the time step, from `dense_overflow_timestep`
- `stop_time`: 40 hours, by which the plume has crossed the slope and spread over the basin floor
- `save_interval`: output interval

A plug of 10 °C water 20 km wide is released on a 500 m shelf into a 20 °C basin 2000 m deep and descends the
slope, held back by the quadratic bottom drag and the lateral viscosity of the benchmark. The configuration
carries no tracer diffusivity, so the whole increase of reference potential energy over the run is spurious,
and the effective diapycnal diffusivity read off it measures the discretization alone.

# Returns
- `Simulation` object after running to completion
"""
function dense_overflow(timestepper::Symbol;
                        arch = CPU(),
                        grid = dense_overflow_grid(arch),
                        free_surface = nothing,
                        averaging_kernel = OptimizedAsymmetricAveragingKernel(),
                        barotropic_timestepper = RungeKutta3Scheme(),
                        slow_forcing = FrozenSlowForcing(),
                        free_surface_name = nothing,
                        boundary_scheme = nothing,
                        tracer_boundary_scheme = something(boundary_scheme, default_tracer_boundary_scheme),
                        momentum_boundary_scheme = something(boundary_scheme, default_momentum_boundary_scheme),
                        tracer_advection = nothing,
                        momentum_advection = nothing,
                        Δt = dense_overflow_timestep(Val(timestepper)),
                        stop_time = 40hours,
                        save_interval = 30minutes)

    p = dense_overflow_parameters()

    tracer_boundary_scheme   = boundary_scheme_value(tracer_boundary_scheme)
    momentum_boundary_scheme = boundary_scheme_value(momentum_boundary_scheme)

    momentum_advection = something(momentum_advection,
                                   split_momentum_advection(nothing, ExplicitTimeDiscretization(),
                                                            momentum_boundary_scheme))

    tracer_advection = something(tracer_advection,
                                 tracer_advection_scheme(tracer_boundary_scheme, TimestepperTestCases.tracer_advection))

    if isnothing(free_surface)
        substeps = dense_overflow_substeps(barotropic_timestepper, timestepper; averaging_kernel)
        free_surface = SplitExplicitFreeSurface(grid; substeps, averaging_kernel,
                                                timestepper = barotropic_timestepper, slow_forcing)
    end

    equation_of_state = LinearEquationOfState(thermal_expansion = p.α, haline_contraction = p.β)

    u_drag = FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form = true, parameters = p.Cᴰ)
    v_drag = FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form = true, parameters = p.Cᴰ)
    u_immersed = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form = true, parameters = p.Cᴰ)
    v_immersed = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form = true, parameters = p.Cᴰ)

    u_bcs = FieldBoundaryConditions(bottom = u_drag, immersed = ImmersedBoundaryCondition(bottom = u_immersed))
    v_bcs = FieldBoundaryConditions(bottom = v_drag, immersed = ImmersedBoundaryCondition(bottom = v_immersed))

    model = HydrostaticFreeSurfaceModel(grid;
                                        timestepper,
                                        coriolis = nothing,
                                        tracers = (:T, :S),
                                        buoyancy = SeawaterBuoyancy(; equation_of_state),
                                        closure = HorizontalScalarDiffusivity(ν = p.ν, κ = 0),
                                        boundary_conditions = (u = u_bcs, v = v_bcs),
                                        free_surface,
                                        momentum_advection,
                                        tracer_advection)

    Tᵢ(x, z) = ifelse(x < p.Lᵖ, p.Tᵖ, p.Tᵃ)

    set!(model, T = Tᵢ, S = p.S)

    simulation = Simulation(model; Δt, stop_time)

    add_callback!(simulation, print_progress, IterationInterval(200))

    ϵb = TimestepperTestCases.BuoyancyVarianceDissipationComputations.BuoyancyVarianceDissipation(grid)
    fb = TimestepperTestCases.BuoyancyVarianceDissipationComputations.flatten_dissipation_fields(ϵb)
    add_callback!(simulation, ϵb, IterationInterval(1))

    VFCC = Oceananigans.AbstractOperations.grid_metric_operation((Face,   Center, Center), Oceananigans.Operators.volume, grid)
    VCCF = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Face),   Oceananigans.Operators.volume, grid)
    VCCC = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Center), Oceananigans.Operators.volume, grid)

    b = Oceananigans.Models.buoyancy_operation(model)
    u, v, w = model.velocities
    T, S = model.tracers
    η = model.free_surface.displacement

    outputs = merge((; u = u * VFCC,
                       w = w * VCCF,
                       T = T * VCCC,
                       S = S * VCCC,
                       b = b * VCCC,
                       Gbx = ∂x(b)^2 * VFCC,
                       Gbz = ∂z(b)^2 * VCCF), fb)

    fsname = something(free_surface_name, default_free_surface_name(free_surface)) *
             boundary_scheme_suffix(tracer_boundary_scheme, momentum_boundary_scheme)

    simulation.output_writers[:fields] = JLD2Writer(model, merge(outputs, (; η));
                                                    filename = "dense_overflow_$(fsname)_$(string(timestepper))",
                                                    schedule = TimeInterval(save_interval),
                                                    overwrite_existing = true)

    run!(simulation)

    return simulation
end

"""
    dense_overflow(d::Discretization; arch = CPU(), kw...)

Run the dense overflow with the discretization `d` of Table 1.
"""
function dense_overflow(d::Discretization; arch = CPU(), kw...)
    return dense_overflow(d.timestepper; arch,
                          barotropic_timestepper = d.barotropic_timestepper,
                          averaging_kernel = d.averaging_kernel,
                          slow_forcing = d.slow_forcing,
                          free_surface = d.implicit_free_surface ? ImplicitFreeSurface() : nothing,
                          tracer_advection = forwarded_tracer_advection(d),
                          free_surface_name = d.label, kw...)
end
