using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids
using Oceananigans.ImmersedBoundaries: immersed_inactive_node
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity
using Random

"""
    wind_stress(i, j, grid, clock, fields, p)

Compute the wind stress forcing at grid point `(i, j)`.

$(SIGNATURES)

# Arguments
- `i, j`: Grid indices
- `grid`: Grid object
- `clock`: Simulation clock
- `fields`: Model fields
- `p`: Parameters named tuple containing `τ₀` (stress amplitude) and `f` (frequency)

# Returns
- Wind stress [m²/s²] that activates after 4 days: `τ₀ * sin(f * t)` if `t > 4 days`, else `0`

This function implements a time-dependent wind stress that begins after a 4-day spin-up period
and oscillates with frequency `f` and amplitude `τ₀`.
"""
@inline function wind_stress(i, j, grid, clock, fields, p) 
    force = clock.time > 4days
    τx = p.τ₀ * sin(p.f * clock.time)
    return ifelse(force, τx, zero(grid))
end

"""
    idealized_coast_timestep(::Val{timestepper})

Return the recommended time step for the idealized coast test case given a timestepper.

$(SIGNATURES)

# Arguments
- `timestepper`: Symbol indicating the timestepper (`:QuasiAdamsBashforth2` or `:SplitRungeKutta3`)

# Returns
- Recommended time step [s] for the given timestepper

Time steps are chosen to match computational cost between AB2 and RK schemes.
"""
idealized_coast_timestep(::Val{scheme}; lowres = false) where scheme =
    baroclinic_timestep(scheme; idealized_coast_stability_parameters(; lowres)...)

"""
    idealized_coast_stability_parameters(; lowres = false)

Stratification, depth and grid spacing that set the first baroclinic Courant number of this case. The shelf is
only 103 m deep, so `c₁` is a fraction of a metre per second and the baroclinic limit is loose; the barotropic
sub-cycle, not the outer scheme, is what the substep count has to keep in hand.
"""
function idealized_coast_stability_parameters(; lowres = false)
    Lx = 192kilometers
    Lz = 103meters
    Nx = lowres ? 96 : 250
    # The shelf is shallow, so c₁ is a fraction of a metre per second and the wave limit is loose; the eddying
    # flow, not the first baroclinic wave, is what sets the time step here. `U` is the horizontal speed the
    # adjustment reaches, and at 1.5 m/s the advective limit is about 4.6 times tighter than the wave one.
    return (N² = 1e-4, H = Lz, Δx = Lx / Nx, U = 1.5, safety = 1)
end

"""
    idealized_coast_substeps(barotropic_scheme, scheme; lowres, averaging_kernel)

Barotropic substep count for this case, fixed so that the substep Courant number sits at 70% of the limit of
`barotropic_scheme` -- `√3` for the three-stage Runge-Kutta substep, `1` for forward-backward.
"""
function idealized_coast_substeps(barotropic_scheme, scheme = :SplitRungeKutta3;
                                  lowres = false, averaging_kernel = OptimizedAsymmetricAveragingKernel())
    p  = idealized_coast_stability_parameters(; lowres)
    Δt = idealized_coast_timestep(Val(scheme); lowres)
    return barotropic_substeps(barotropic_scheme; p.H, p.Δx, Δt, averaging_kernel)
end

@inline ϕ²(i, j, k, grid, ϕ)    = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_quadratic_bottom_drag(i, j, grid, clock, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_quadratic_bottom_drag(i, j, grid, clock, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, fields)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, fields)

"""
    idealized_coast(timestepper::Symbol; arch, forced, lowres, free_surface)

Set up and run the idealized coastal baroclinic adjustment test case simulation.

$(SIGNATURES)

# Arguments
- `timestepper`: Symbol indicating the timestepper (`:QuasiAdamsBashforth2` or `:SplitRungeKutta3`)

# Keyword Arguments
- `arch`: Architecture to run on (default: `CPU()`)
- `forced`: Whether to apply wind forcing (default: `true`)
- `lowres`: Whether to use low resolution (default: `false`)
- `free_surface`: Free surface formulation (default: `SplitExplicitFreeSurface` with CFL=0.7)

# Returns
- `Simulation` object configured but not yet run

This function sets up the idealized coastal baroclinic adjustment test case described in the paper.
The configuration features a rectangular domain with a linearly sloping bathymetry and initial
meridional salinity gradient that drives baroclinic instabilities. The test case demonstrates
how numerical mixing can suppress submesoscale variability, as shown in the paper.

The simulation runs for 40 days and outputs velocity, temperature, salinity, buoyancy,
and variance dissipation diagnostics.
"""

function idealized_coast(timestepper::Symbol;
                         arch = CPU(),
                         forced = false,
                         lowres = false,
                         free_surface = nothing,
                         averaging_kernel = OptimizedAsymmetricAveragingKernel(),
                         barotropic_timestepper = ForwardBackwardScheme(),
                         slow_forcing = FrozenSlowForcing(),
                         free_surface_name = nothing,
                         tracer_advection = TimestepperTestCases.tracer_advection,
                         stop_time = 40days,
                         stop_iteration = Inf)

    Lx = 192kilometers
    Ly = 192kilometers
    Lz = 103meters

    if lowres
        Nx = Ny = 96
    else
        Nx = Ny = 250
    end
        
    Nz = 60   # 103 m over 60 levels is the 1.7 m vertical resolution the manuscript quotes
    z_faces = (-Lz, 0)
    
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
    Δt = idealized_coast_timestep(Val(timestepper); lowres)

    if isnothing(free_surface)
        # Fixed substep count, set from the substep integrator's own stability limit rather than from a cfl
        # target, so that the barotropic Courant number is a stated property of the run.
        substeps = idealized_coast_substeps(barotropic_timestepper, timestepper; lowres, averaging_kernel)
        free_surface = SplitExplicitFreeSurface(grid; substeps, averaging_kernel,
                                                timestepper=barotropic_timestepper, slow_forcing)
    end

    τ₀ = 0.1 / 1027

    bottom_drag_coefficient = 0.003

    u_immersed_drag = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
    v_immersed_drag = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)

    u_immersed = ImmersedBoundaryCondition(bottom=u_immersed_drag)
    v_immersed = ImmersedBoundaryCondition(bottom=v_immersed_drag)
    u_bottom   = FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
    v_bottom   = FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
    u_top      = FluxBoundaryCondition(wind_stress; discrete_form=true, parameters=(τ₀=τ₀, f=f))
   
    u_bcs = if forced
        FieldBoundaryConditions(bottom=u_bottom, immersed=u_immersed, top=u_top)
    else
        FieldBoundaryConditions(bottom=u_bottom, immersed=u_immersed)
    end
    v_bcs = FieldBoundaryConditions(bottom=v_bottom, immersed=v_immersed)

    cl1 = forced ? ConvectiveAdjustmentVerticalDiffusivity(background_κz=1e-5, 
                                                           convective_κz=0.1, 
                                                           background_νz=1e-5, 
                                                           convective_νz=0.1) : nothing

    tracers = (:T, :S) 

    model = HydrostaticFreeSurfaceModel(grid;
                                        coriolis,
                                        timestepper,
                                        tracers,
                                        buoyancy,
                                        closure = cl1,
                                        boundary_conditions = (; u=u_bcs, v=v_bcs),
                                        free_surface,
                                        momentum_advection = WENOVectorInvariant(),
                                        tracer_advection)

    N² = 1e-4
    S² = 1e-8
    M²(y) = if y > 50kilometers
        0.0
    else
        1.2e-6
    end
    g  = buoyancy.gravitational_acceleration

    Random.seed!(1234)

    Tᵢ(x, y, z) = 25 + N² / (α * g) * z + 1e-2 * rand()
    Sᵢ(x, y, z) = 35 - M²(y) / (β * g) * (50kilometers - y) - S² / (β * g) * z
    uᵢ(x, y, z) = y > 60kilometers ? 0.0 : - 1 / f * M²(y) * (z - bottom_height(x, y))

    set!(model, T=Tᵢ, S=Sᵢ)
    simulation = Simulation(model; Δt, stop_time, stop_iteration)

    add_callback!(simulation, print_progress,  IterationInterval(100))

    # Dissipations...
    ϵT = Oceananigans.Models.VarianceDissipationComputations.VarianceDissipation(:T, grid)
    ϵS = Oceananigans.Models.VarianceDissipationComputations.VarianceDissipation(:S, grid)
    fT = Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵT)
    fS = Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵS)
    add_callback!(simulation, ϵT, IterationInterval(1))
    add_callback!(simulation, ϵS, IterationInterval(1))
    
    ϵb = TimestepperTestCases.BuoyancyVarianceDissipationComputations.BuoyancyVarianceDissipation(grid)
    fb = TimestepperTestCases.BuoyancyVarianceDissipationComputations.flatten_dissipation_fields(ϵb)
    add_callback!(simulation, ϵb, IterationInterval(1))

    if free_surface_name !== nothing
        fsname = free_surface_name
    elseif free_surface isa SplitExplicitFreeSurface
        fsname = "split_free_surface"
    else
        fsname = "implicit_free_surface"
    end

    filename = "idealized_coast_$(fsname)_$(lowres ? "lowres" : "")"
    save_fields_interval = 12hours

    closure = cl1 isa CATKEVerticalDiffusivity ? "CATKE" : cl1 isa Nothing ? "unforced" : "RiBased"
   
    VFCC = Oceananigans.AbstractOperations.grid_metric_operation((Face,   Center, Center), Oceananigans.Operators.volume, grid)
    VCFC = Oceananigans.AbstractOperations.grid_metric_operation((Center, Face,   Center), Oceananigans.Operators.volume, grid)
    VCCF = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Face),   Oceananigans.Operators.volume, grid)
    VCCC = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Center), Oceananigans.Operators.volume, grid)

    T, S = model.tracers

    # GTx = ∂x(T)^2 * VFCC
    # GTy = ∂y(T)^2 * VCFC
    # GTz = ∂z(T)^2 * VCCF
    # GSx = ∂x(S)^2 * VFCC
    # GSy = ∂y(S)^2 * VCFC
    # GSz = ∂z(S)^2 * VCCF

    b   = Oceananigans.Models.buoyancy_operation(model)
    Gbx = ∂x(b)^2 * VFCC
    Gby = ∂y(b)^2 * VCFC
    Gbz = ∂z(b)^2 * VCCF

    G = (; Gbx, Gby, Gbz) # GTx, GTy, GTz, GSx, GSy, GSz)
    u, v, w = model.velocities
    η = model.free_surface.displacement
    T, S = model.tracers

    outputs = merge((; u = u * VFCC,
                       v = v * VCFC,
                       w = w * VCCF,
                       T = T * VCCC,
                       S = S * VCCC,
                       b = b * VCCC), fb, G) # , fT, fS

    if !isnothing(cl1)
        κu = model.closure_fields.κu
        κc = model.closure_fields.κc
        outputs = merge(outputs, (; κu, κc))
    end

    average_outputs = NamedTuple{keys(outputs)}(Average(output, dims=1) for output in values(outputs))

    simulation.output_writers[:values] = JLD2Writer(model, merge(outputs, (; η));
                                                    filename = filename * "_$(string(timestepper))_$(closure)",
                                                    schedule = TimeInterval(save_fields_interval),
                                                    file_splitting = TimeInterval(1days),
                                                    array_type = Array{Float32},
                                                    overwrite_existing = true)

    simulation.output_writers[:surface] = JLD2Writer(model, outputs;
                                                    filename = filename * "_surface_$(string(timestepper))_$(closure)",
                                                    schedule = TimeInterval(save_fields_interval),
                                                    array_type = Array{Float32},
                                                    indices = (:, :, grid.Nz),
                                                    overwrite_existing = true)

    simulation.output_writers[:averages] = JLD2Writer(model, average_outputs;
                                                      filename = filename * "_averages_$(string(timestepper))_$(closure)",
                                                      schedule = TimeInterval(save_fields_interval),
                                                      array_type = Array{Float32},
                                                      overwrite_existing = true)

    @info "Running the simulation..."

    run!(simulation)

    return simulation
end

"""
    idealized_coast(d::Discretization; arch = CPU(), lowres = false, kw...)

Run the coastal adjustment case with the discretization `d` of Table 1.
"""
function idealized_coast(d::Discretization; arch = CPU(), lowres = false, kw...)
    return idealized_coast(d.timestepper; arch, lowres,
                           barotropic_timestepper = d.barotropic_timestepper,
                           averaging_kernel = d.averaging_kernel,
                           slow_forcing = d.slow_forcing,
                           free_surface = d.implicit_free_surface ? ImplicitFreeSurface() : nothing,
                           tracer_advection = d.tracer_advection,
                           free_surface_name = d.label, kw...)
end
