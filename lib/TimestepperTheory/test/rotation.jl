function rotation_checks()
    modes = vertical_modes(0.05, 8)
    N = length(modes.α)
    indices = [collect(1:N); collect(2N+1:3N)]
    filter = averaging_filter(F3Filter(), 96)
    for (coupling, scheme, forcing) in ((FrozenRotation(), WickerSkamarock3(), FrozenForcing()),
                                       (ProgressivePressure(), WickerSkamarock3(), ProgressiveForcing()),
                                       (ProgressiveRotation(), WickerSkamarock3(), ProgressiveForcing()),
                                       (SSPRotation(), ShuOsher3(), FrozenForcing()))
        M = rotating_evolution_matrix(2., 0., modes, coupling; filter)
        M₀ = evolution_matrix(2., modes; scheme, forcing, filter)
        @test M[indices, indices] ≈ M₀ atol=2e-13
        @test M[N+1:2N, N+1:2N] ≈ I atol=2e-13
        @test sort(abs.(eigvals(rotating_evolution_matrix(2., 0.3, modes, coupling; filter)))) ≈
              sort(abs.(eigvals(rotating_evolution_matrix(2., -0.3, modes, coupling; filter)))) atol=2e-12
    end
end


@testset "Rotating modal equations" begin
    rotation_checks()
end

@testset "SRK3 inertial amplification" begin
    filter = averaging_filter(F3Filter(), 96)
    modes = (; α=[1.], χ=[1.], ψ=[1.], Π=[0.])
    for ρ in (0.1, 0.3, 0.5)
        M = rotating_evolution_matrix(1e-7, ρ, modes, SSPRotation(); filter)
        λ = 1 + im*ρ - ρ^2/2 - im*ρ^3/6
        @test M[1:2, 1:2] ≈ [real(λ) imag(λ); -imag(λ) real(λ)] atol=2e-13
    end
end

@testset "Rotating manuscript configurations" begin
    modes = vertical_modes(0.05, 8)
    N = length(modes.α)
    indices = [collect(1:N); collect(2N+1:3N)]
    for (scheme, coupling, forcing, averaging) in ((ModifiedRungeKutta4(a=0.23), ProgressiveRotation(), ProgressiveForcing(), F3Filter()),
                                                  (WickerSkamarock3(), FrozenRotation(), FrozenForcing(), F1Filter()))
        filter = averaging_filter(averaging, 96)
        M = rotating_evolution_matrix(2., 0., modes, coupling; scheme, filter)
        @test M[indices, indices] ≈ evolution_matrix(2., modes; scheme, forcing, filter) atol=2e-13
        @test M[N+1:2N, N+1:2N] ≈ I atol=2e-13
    end
end
