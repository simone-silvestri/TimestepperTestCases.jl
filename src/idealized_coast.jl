using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids
using Oceananigans.ImmersedBoundaries: immersed_inactive_node
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity

@inline function wind_stress(i, j, grid, clock, fields, p) 
    force = clock.time > 3days
    τx = p.τ₀ * sin(0.92 * p.f * clock.time)
    return ifelse(force, τx, zero(grid))
end

idealized_coast_timestep(::Val{:QuasiAdamsBashforth2}) = 60seconds
idealized_coast_timestep(::Val{:SplitRungeKutta3})     = 60seconds #3minutes

@inline ϕ²(i, j, k, grid, ϕ)    = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, fields)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, fields)

function idealized_coast(timestepper::Symbol; arch = CPU(), forced = true)

    Lx = 97kilometers
    Ly = 97kilometers
    Lz = 103meters

    Nx = 194
    Ny = 194
    Nz = 120

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

    bottom_drag_coefficient = 0.03

    u_immersed_drag = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
    v_immersed_drag = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)

    u_immersed = ImmersedBoundaryCondition(bottom=u_immersed_drag)
    v_immersed = ImmersedBoundaryCondition(bottom=v_immersed_drag)
    u_bottom   = FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
    v_bottom   = FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
    u_top      = FluxBoundaryCondition(wind_stress; discrete_form=true, parameters=(τ₀=τ₀, f=f))
   
    u_bcs = FieldBoundaryConditions(bottom=u_bottom, immersed=u_immersed, top=u_top)
    v_bcs = FieldBoundaryConditions(bottom=v_bottom, immersed=v_immersed)

    # cl1 = ConvectiveAdjustmentVerticalDiffusivity(convective_κz=0.1, convective_νz=0.1) 
    cl1 = RiBasedVerticalDiffusivity(horizontal_Ri_filter=Oceananigans.TurbulenceClosures.FivePointHorizontalFilter())
    cl2 = VerticalScalarDiffusivity(ν=1e-4)

    model = HydrostaticFreeSurfaceModel(; grid,
                                          coriolis,
                                          timestepper,
                                          tracers = (:T, :S),
                                          buoyancy,
                                          closure = (cl1, cl2),
                                          boundary_conditions = (; u=u_bcs, v=v_bcs),
                                          free_surface = SplitExplicitFreeSurface(grid; substeps=50),
                                          momentum_advection = WENOVectorInvariant(),
                                          tracer_advection = WENO(; order = 7, buffer_scheme = Centered()))

    N² = 1e-4
    

    N² = 1e-4
    M²(y) = if y > 50kilometers
        0.0
    else
        1e-6
    end
    g  = buoyancy.gravitational_acceleration

    Tᵢ(x, y, z) = 25 + N² / (α * g) * z
    Sᵢ(x, y, z) = 35 - M²(y) / (β * g) * (50kilometers - y)
    uᵢ(x, y, z) = y > 60kilometers ? 0.0 : - 1 / f * M²(y) * (z - bottom_height(x, y))

    set!(model, T=Tᵢ, S=Sᵢ, u=uᵢ)
    Δt = idealized_coast_timestep(Val(timestepper))
    simulation = Simulation(model; Δt, stop_time=30days)

    add_callback!(simulation, print_progress,  IterationInterval(100))

    filename = "idealized_coast"
    save_fields_interval = 1hours

    simulation.output_writers[:values] = JLD2Writer(model, merge(model.velocities, model.tracers);
                                                    filename = filename * "_$(string(timestepper))",
                                                    schedule = IterationInterval(300),
                                                    overwrite_existing = true)


    @info "Running the simulation..."

    if timestepper == :SplitRungeKutta3
        conjure_time_step_wizard!(simulation, cfl=0.7, IterationInterval(10), max_Δt=3minutes)
    end
    
#    run!(simulation)

    @info "Simulation completed in " * prettytime(simulation.run_wall_time)

    return simulation
end

