using TimestepperTestCases
using Oceananigans
using Oceananigans.Grids: topology
using Oceananigans.Units
using Test

@testset "configuration" begin
    p = dense_overflow_parameters()

    # the shelf and the basin are the two plateaux of the tanh, reached well outside the slope
    @test TimestepperTestCases.dense_overflow_bottom_height(0.0) ≈ -p.Hˢ atol = 1
    @test TimestepperTestCases.dense_overflow_bottom_height(p.Lx) ≈ -p.H atol = 1
    @test TimestepperTestCases.dense_overflow_bottom_height(p.xˢ) ≈ -(p.H + p.Hˢ) / 2

    grid = dense_overflow_grid()
    @test size(grid) == (p.Nx, 1, p.Nz)
    @test topology(grid) == (Bounded, Flat, Bounded)

    # every scheme sits at the same fraction of its own limit, so the steps are in the ratio of the limits
    Δt₃ = dense_overflow_timestep(Val(:SplitRungeKutta3))
    Δt₂ = dense_overflow_timestep(Val(:QuasiAdamsBashforth2))
    @test Δt₂ / Δt₃ ≈ stability_limit(:QuasiAdamsBashforth2) / stability_limit(:SplitRungeKutta3)
end

@testset "integration" begin
    simulation = dense_overflow(:SplitRungeKutta3; stop_time = 1hour, save_interval = 30minutes)
    model = simulation.model

    T = interior(model.tracers.T)
    S = interior(model.tracers.S)
    wet = T .!= 0

    p = dense_overflow_parameters()

    # the plug is bounded by the two water masses it was released between
    @test minimum(T[wet]) ≥ p.Tᵖ - 0.1
    @test maximum(T[wet]) ≤ p.Tᵃ + 0.1

    # salinity is uniform and carries no diffusivity, so advecting it must leave it alone
    @test all(≈(p.S), S[wet])

    @test all(isfinite, interior(model.velocities.u))
end

@testset "reference potential energy" begin
    case = load_dense_overflow(".", "split_free_surface", "SplitRungeKutta3")

    # mixing is irreversible and the case carries no diffusivity, so the whole rise is spurious
    @test issorted(case[:rpe])
    @test dense_overflow_diffusivity(case) > 0
end
