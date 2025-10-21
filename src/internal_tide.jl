using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids

function internal_tide_parameters() 
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

@inline tidal_forcing(x, z, t, p) = p.A₂ * sin(p.ω₂ * t)

internal_tide_timestep(::Val{:QuasiAdamsBashforth2}) = 5minutes
internal_tide_timestep(::Val{:SplitRungeKutta3})     = 15minutes
internal_tide_timestep(::Val{:SplitRungeKutta4})     = 20minutes

@kernel function _compute_dissipation!(Δtc², c⁻, c, Δt)
    i, j, k = @index(Global, NTuple)
    @inbounds begin
        Δtc²[i, j, k] = (c[i, j, k]^2 - c⁻[i, j, k]^2) / Δt
        c⁻[i, j, k]   = c[i, j, k]
    end
end

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

function internal_tide_grid()
    param = internal_tide_parameters()

    Nx, Nz    = param.Nx, param.Nz
    h₀, width = param.h₀, param.width
    H, L      = param.H, param.L

    underlying_grid = RectilinearGrid(size = (Nx, Nz), halo = (6, 6),
                                    x = (-L, L), z = MutableVerticalDiscretization((-H, 0)),
                                    topology = (Periodic, Flat, Bounded))

    hill(x)   =   h₀ * exp(-x^2 / 2width^2)
    bottom(x) = - H + hill(x)

    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom); active_cells_map = true)
   
    return grid
end

function internal_tide(timestepper::Symbol; free_surface=SplitExplicitFreeSurface(internal_tide_grid(); substeps=60))
    
    grid  = internal_tide_grid()
    param = internal_tide_parameters()

    coriolis  = FPlane(f = param.f)
    u_forcing = Forcing(tidal_forcing, parameters=param)

    c⁻    = CenterField(grid)
    Δtc²  = CenterField(grid)
    b⁻    = CenterField(grid)
    Δtb²  = CenterField(grid)

    model = HydrostaticFreeSurfaceModel(; grid, coriolis,
                                          buoyancy = BuoyancyTracer(),
                                          tracers = (:b, :c),
                                          momentum_advection = WENO(),
                                          tracer_advection,
                                          free_surface,
                                          timestepper,
                                          forcing = (; u = u_forcing),
                                          auxiliary_fields=(; Δtc², c⁻, Δtb², b⁻))


    bᵢ(x, z) = param.Nᵢ² * z
    cᵢ(x, z) = exp( - (z + 1kilometers)^2 / (2 * (25meters)^2))
    set!(model, u=param.U₂, b=bᵢ, c=cᵢ)

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
        fsname = "split_free_surface"
    else
        fsname = "implicit_free_surface"
    end

    filename = "internal_tide_$(string(timestepper))_$(round(Δt/minutes))min_$(fsname)"
    save_fields_interval = 1hours
    
    f = merge(Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵb),
              Oceananigans.Models.VarianceDissipationComputations.flatten_dissipation_fields(ϵc))

    f = (; Abx = f.Abx,
           Abz = f.Abz,
           Acx = f.Acx,
           Acz = f.Acz)

    VFC = Oceananigans.AbstractOperations.grid_metric_operation((Face,   Center, Center), Oceananigans.Operators.volume, grid)
    VCF = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Face),   Oceananigans.Operators.volume, grid)
    VCC = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Center), Oceananigans.Operators.volume, grid)

    V  = (; VFC, VCF, VCC)
    Δ² = (; Δtc² = model.auxiliary_fields.Δtc², Δtb² = model.auxiliary_fields.Δtb²)

    simulation.output_writers[:fields] = JLD2Writer(model, merge((; u, u′, w, b, c, N², η), f, g, V, Δ²); filename,
                                                    schedule = TimeInterval(save_fields_interval),
                                                    overwrite_existing = true)

    run!(simulation)

    return simulation
end

function visualize_internal_tide(filename)
    u′_t = FieldTimeSeries(filename, "u′")
     w_t = FieldTimeSeries(filename, "w")
     c_t = FieldTimeSeries(filename, "c")
    N²_t = FieldTimeSeries(filename, "N²")
    ϵx_t = FieldTimeSeries(filename, "Abx")
    ϵz_t = FieldTimeSeries(filename, "Abz")

    umax = maximum(abs, u′_t[end])
    wmax = maximum(abs, w_t[end])

    grid  = u′_t.grid
    times = u′_t.times

    n = Observable(1)
    param = internal_tide_parameters()
    title = @lift @sprintf("t = %1.2f days = %1.2f T₂",
                        round(times[$n] / day, digits=2) , round(times[$n] / param.T₂, digits=2))

    u′ₙ = @lift interior(u′_t[$n], :, 1, :)
    wₙ =  @lift interior( w_t[$n], :, 1, :)
    cₙ =  @lift interior( c_t[$n], :, 1, :)
    N²ₙ = @lift interior(N²_t[$n], :, 1, :)
    ϵxₙ = @lift interior(ϵx_t[$n], :, 1, :)
    ϵzₙ = @lift interior(ϵz_t[$n], :, 1, :)

    axis_kwargs = (xlabel = "x [m]",
                   ylabel = "z [m]",
                   limits = ((-grid.Lx/2, grid.Lx/2), (-grid.Lz, 0)),
                   titlesize = 20)

    fig = Figure(size = (700, 1500))

    fig[1, :] = Label(fig, title, fontsize=24, tellwidth=false)

    ax_u = Axis(fig[2, 1]; title = "u'-velocity") #, axis_kwargs...)
    hm_u = heatmap!(ax_u, u′ₙ; nan_color=:gray, colorrange=(-0.35, 0.35), colormap=:balance)
    Colorbar(fig[2, 2], hm_u, label = "m s⁻¹")

    ax_w = Axis(fig[3, 1]; title = "w-velocity") #, axis_kwargs...)
    hm_w = heatmap!(ax_w, wₙ; nan_color=:gray, colorrange=(-0.0035, 0.0035), colormap=:balance)
    Colorbar(fig[3, 2], hm_w, label = "m s⁻¹")

    ax_N² = Axis(fig[4, 1]; title = "stratification N²")# , axis_kwargs...)
    hm_N² = heatmap!(ax_N², N²ₙ; nan_color=:gray, colorrange=(0.9param.Nᵢ², 1.1param.Nᵢ²), colormap=:magma)
    Colorbar(fig[4, 2], hm_N², label = "s⁻²")

    ax_ϵx = Axis(fig[5, 1]; title = "variance dissipation rate ϵ") #, axis_kwargs...)
    hm_ϵx = heatmap!(ax_ϵx, ϵxₙ; nan_color=:gray, colorrange=(-3e-7, 3e-7), colormap=:magma)
    Colorbar(fig[5, 2], hm_ϵx, label = "m² s⁻³")

    ax_ϵz = Axis(fig[6, 1]; title = "variance dissipation rate ϵ_z") #, axis_kwargs...)
    hm_ϵz = heatmap!(ax_ϵz, ϵzₙ; nan_color=:gray, colorrange=(-3e-7, 3e-7), colormap=:magma)
    Colorbar(fig[6, 2], hm_ϵz, label = "m² s⁻³")

    ax_c = Axis(fig[7, 1]; title = "tracer c") #, axis_kwargs...)
    hm_c = heatmap!(ax_c, cₙ; nan_color=:gray, colorrange=(0, 1), colormap=:viridis)
    Colorbar(fig[7, 2], hm_c, label = "concentration")

    fig

    @info "Making an animation from saved data..."

    frames = 1:length(times)

    record(fig, filename * ".mp4", frames, framerate=16) do i
        @info string("Plotting frame ", i, " of ", frames[end])
        n[] = i
    end

    return fig
end

function compute_kinetic_energy(filename)
    u   = FieldTimeSeries(filename, "u")
    w   = FieldTimeSeries(filename, "w")
    VFC = FieldTimeSeries(filename, "VFC") 
    VCF = FieldTimeSeries(filename, "VCF")
    KE = zeros(length(u))

    for i in 1:length(u)
        ut = u[i]
        wt = w[i]
        KE[i] = 0.5 * mean(ut^2 * VFC[i] + wt^2 * VCF[i])
    end

    return KE
end

function compute_budget_dissipation(filename, var)
    Ax  = FieldTimeSeries(filename, "A" * var * "x")
    Az  = FieldTimeSeries(filename, "A" * var * "z")
    Δ²  = FieldTimeSeries(filename, "Δ" * var * "²")
    
    VFC = FieldTimeSeries(filename, "VFC")
    VCF = FieldTimeSeries(filename, "VCF")
    VCC = FieldTimeSeries(filename, "VCC")
    
    ∫Ax = zeros(length(Ax))
    ∫Az = zeros(length(Az))
    ∫Δ² = zeros(length(Δ²))
    
    for i in 1:length(Ax)
        ∫Ax[i] = abs(sum(Az[i] * VFC[i]))
        ∫Az[i] = abs(sum(Az[i] * VCF[i]))
        ∫Δ²[i] = abs(sum(Δ²[i] * VCC[i]))
    end

    return (; ∫Ax, ∫Az, ∫Δ²)
end

function compute_dissipation(filename)
    b   = FieldTimeSeries(filename, "b")
    b2  = zeros(length(b))
    VCC = FieldTimeSeries(filename, "VCC")
    for i in 1:length(b)
        b2[i] = mean(b[i]^2 * VCC[i])
    end
    return b2
end

function compute_conservation(filename, var)
    b   = FieldTimeSeries(filename, var)
    bi  = zeros(length(b))
    VCC = FieldTimeSeries(filename, "VCC")

    for i in 1:length(b)
        bi[i] = sum(b[i] * VCC[i]) / sum(VCC[i])
    end
    
    return bi
end
