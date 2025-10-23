using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids
using Oceananigans.ImmersedBoundaries: immersed_inactive_node
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity

@inline function wind_stress(i, j, grid, clock, fields, p) 
    force = clock.time > 4days
    τx = p.τ₀ * sin(p.f * clock.time)
    return ifelse(force, τx, zero(grid))
end

idealized_coast_timestep(::Val{:QuasiAdamsBashforth2}) = 3minutes
idealized_coast_timestep(::Val{:SplitRungeKutta3})     = 9minutes

@inline ϕ²(i, j, k, grid, ϕ)    = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, fields)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, fields)

function idealized_coast(timestepper::Symbol; 
                         arch = CPU(), 
                         forced = true,
                         closure = CATKEVerticalDiffusivity(),
                         free_surface = nothing)

    Lx = 96kilometers
    Ly = 96kilometers
    Lz = 103meters

    Nx = 250
    Ny = 250
    Nz = 60

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

    if isnothing(free_surface)
        free_surface = SplitExplicitFreeSurface(grid; substeps=50)
    end

    coriolis = FPlane(; f)

    equation_of_state = LinearEquationOfState(thermal_expansion=α, haline_contraction=β)
    buoyancy = SeawaterBuoyancy(; equation_of_state)

    if isnothing(free_surface)
        free_surface = SplitExplicitFreeSurface(grid; substeps=60)
    end

    if forced 
        τ₀ = 0.05 / 1027
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

    cl1 = closure 
    cl2 = VerticalScalarDiffusivity(ν=1e-4)

    model = HydrostaticFreeSurfaceModel(; grid,
                                          coriolis,
                                          timestepper,
                                          tracers = (:T, :S, :e),
                                          buoyancy,
                                          closure = (cl1, cl2),
                                          boundary_conditions = (; u=u_bcs, v=v_bcs),
                                          free_surface,
                                          momentum_advection = WENOVectorInvariant(),
                                          tracer_advection)

    N² = 1e-4
    S² = 1e-6
    M²(y) = if y > 50kilometers
        0.0
    else
        1.5e-6
    end
    g  = buoyancy.gravitational_acceleration

    Tᵢ(x, y, z) = 25 + N² / (α * g) * z
    Sᵢ(x, y, z) = 35 - M²(y) / (β * g) * (50kilometers - y) + S² / (β * g) * z
    uᵢ(x, y, z) = y > 60kilometers ? 0.0 : - 1 / f * M²(y) * (z - bottom_height(x, y))

    set!(model, T=Tᵢ, S=Sᵢ, u=uᵢ)
    Δt = idealized_coast_timestep(Val(timestepper))
    simulation = Simulation(model; Δt, stop_time=30days)

    add_callback!(simulation, print_progress,  IterationInterval(100))

    # Dissipations...
    ϵT = Oceananigans.Models.VarianceDissipationComputations.VarianceDissipation(:T, grid)
    ϵS = Oceananigans.Models.VarianceDissipationComputations.VarianceDissipation(:S, grid)
    fT = Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵT)
    fS = Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵS)

    add_callback!(simulation, ϵT, IterationInterval(1))
    add_callback!(simulation, ϵS, IterationInterval(1))

    if free_surface isa SplitExplicitFreeSurface
        fsname = "split_free_surface"
    else
        fsname = "implicit_free_surface"
    end

    filename = "idealized_coast_$(fsname)"
    save_fields_interval = 1hours

    closure = cl1 isa CATKEVerticalDiffusivity ? "CATKE" : "RiBased"

    VFCC = Oceananigans.AbstractOperations.grid_metric_operation((Face,   Center, Center), Oceananigans.Operators.volume, grid)
    VCFC = Oceananigans.AbstractOperations.grid_metric_operation((Center, Face,   Center), Oceananigans.Operators.volume, grid)
    VCCF = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Face),   Oceananigans.Operators.volume, grid)
    VCCC = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Center), Oceananigans.Operators.volume, grid)

    T, S = model.tracers

    GTx = ∂x(T)^2 * VFCC
    GTy = ∂y(T)^2 * VCFC
    GTz = ∂z(T)^2 * VCCF
    GSx = ∂x(S)^2 * VFCC
    GSy = ∂y(S)^2 * VCFC
    GSz = ∂z(S)^2 * VCCF

    b   = Oceananigans.BuoyancyFormulations.buoyancy(model.buoyancy, grid, (; T = T, S = S))
    Gbx = ∂x(b)^2 * VFCC
    Gby = ∂y(b)^2 * VCFC
    Gbz = ∂z(b)^2 * VCCF

    Abx = g * (α * fT.ATx / VFCC - β * fS.ASx / VFCC) * VFCC
    Aby = g * (α * fT.ATy / VCFC - β * fS.ASy / VCFC) * VCFC
    Abz = g * (α * fT.ATz / VCCF - β * fS.ASz / VCCF) * VCCF

    G = (; GTx, GTy, GTz, GSx, GSy, GSz, Gbx, Gby, Gbz)
    u, v, w = model.velocities
    η = model.free_surface.η
    T, S = model.tracers
    κu = model.diffusivity_fields[1].κu
    κc = model.diffusivity_fields[1].κc

    outputs = merge((; u = u * VFCC,
                       v = v * VCFC,
                       w = w * VCCF,
                       T = T * VCCC,
                       S = S * VCCC,
                       b = b * VCCC,
                       κu, κc), fT, fS, (; Abx, Aby, Abz), G)

    average_outputs = NamedTuple{keys(outputs)}(Average(output, dims=1) for output in values(outputs))

    simulation.output_writers[:values] = JLD2Writer(model, merge(outputs, (; η));
                                                    filename = filename * "_$(string(timestepper))_$(closure)",
                                                    schedule = TimeInterval(save_fields_interval),
                                                    file_splitting = TimeInterval(1days),
                                                    array_type = Array{Float32},
                                                    overwrite_existing = true)

    simulation.output_writers[:averages] = JLD2Writer(model, average_outputs;
                                                      filename = filename * "_averages_$(string(timestepper))_$(closure)",
                                                      schedule = TimeInterval(save_fields_interval),
                                                      file_splitting = TimeInterval(1days),
                                                      array_type = Array{Float32},
                                                      overwrite_existing = true)

    @info "Running the simulation..."

    return simulation
end

