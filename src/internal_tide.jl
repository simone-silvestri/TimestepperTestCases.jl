using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids
using Printf

function internal_tide(timestepper::Symbol)
    Nx, Nz = 256, 128
    H, L = 2kilometers, 1000kilometers

    underlying_grid = RectilinearGrid(size = (Nx, Nz), halo = (6, 6),
                                      x = (-L, L), z = MutableVerticalDiscretization((-H, 0)),
                                      topology = (Periodic, Flat, Bounded))

    h₀ = 250meters
    width = 20kilometers
    hill(x) = h₀ * exp(-x^2 / 2width^2)
    bottom(x) = - H + 0 + hill(x)

    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom))
    coriolis = FPlane(latitude = -45)

    T₂ = 12.421hours
    ω₂ = 2π / T₂ # radians/sec
    ϵ = 0.1 # excursion parameter
    U₂ = ϵ * ω₂ * width
    A₂ = U₂ * (ω₂^2 - coriolis.f^2) / ω₂

    @inline tidal_forcing(x, z, t, p) = p.A₂ * sin(p.ω₂ * t)
    u_forcing = Forcing(tidal_forcing, parameters=(; A₂, ω₂))
    free_surface = SplitExplicitFreeSurface(grid; substeps=50)
    model = HydrostaticFreeSurfaceModel(; grid, coriolis,
                                        buoyancy = BuoyancyTracer(),
                                        tracers = :b,
                                        momentum_advection = WENO(),
                                        tracer_advection = WENO(order=7), 
                                        free_surface,
                                        timestepper,
                                        forcing = (; u = u_forcing))

    Nᵢ² = 1e-4  
    bᵢ(x, z) = Nᵢ² * z
    set!(model, u=U₂, b=bᵢ)
    Δt = 15minutes
    stop_time = 40days
    simulation = Simulation(model; Δt, stop_time)

    ϵ = Oceananigans.Models.VarianceDissipationComputations.VarianceDissipation(:b, grid)
    add_callback!(simulation, ϵ, IterationInterval(1))
    add_callback!(simulation, progress, IterationInterval(200))

    b = model.tracers.b
    u, v, w = model.velocities
    U = Field(Average(u))
    u′ = u - U
    N² = ∂z(b)

    filename = "internal_tide_$(string(timestepper))_$(round(Δt/minutes))min"
    save_fields_interval = 30minutes
    f = Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵ)

    simulation.output_writers[:fields] = JLD2Writer(model, merge((; u, u′, w, b, N²), f); filename,
                                                    schedule = TimeInterval(save_fields_interval),
                                                    overwrite_existing = true)

    @info "Running the simulation..."
    run!(simulation)

    @info "Simulation completed in " * prettytime(simulation.run_wall_time)

    return simulation
end