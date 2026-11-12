using TimestepperTheory
using LinearAlgebra
using Test

include("direct_steps.jl")

@testset "Vertical modes" begin
    for number_type in (Float64, BigFloat)
        modes = vertical_modes(number_type(1 // 20), 8)
        @test modes.x .* tan.(modes.x) ≈ fill(number_type(1 // 20), 9)
        @test modes.χ ≈ modes.ψ .* modes.α.^2
        @test all(modes.x .> (0:8) .* number_type(π))
        @test all(modes.x .< ((0:8) .+ 1 // 2) .* number_type(π))
        @test sum(abs2, modes.χ) < 1
        for first in 1:9, second in 1:9
            left, right = modes.x[first], modes.x[second]
            integral = (sinc((left - right) / number_type(π)) + sinc((left + right) / number_type(π))) / 2
            @test modes.A[first] * modes.A[second] * integral ≈ (first == second ? 1 : 0) atol=100eps(number_type)
        end
    end
    coarse = vertical_modes(big"0.01", 8)
    fine = vertical_modes(big"0.01", 16)
    @test 6 < (1 - sum(abs2, coarse.χ)) / (1 - sum(abs2, fine.χ)) < 9
    modes = vertical_modes(big"1e-5", 2)
    @test isapprox(modes.χ[1], 1 - big"1e-10" / 90; atol=big"1e-17")
end

@testset "Filter moments and barotropic continuity" begin
    for name in (NoFilter(), F1Filter(), F2Filter(), F3Filter())
        filter = averaging_filter(name; float_type=BigFloat)
        @test sum(filter.w) ≈ 1
        @test sum(weight * index * filter.Δτ for (index, weight) in enumerate(filter.w)) ≈ 1
        if name in (F2Filter(), F3Filter())
            @test abs(sum(weight * (index * filter.Δτ - 1)^2 for (index, weight) in enumerate(filter.w))) < big"1e-65"
        end
        if name == F3Filter()
            @test abs(sum(weight * (index * filter.Δτ - 1)^3 for (index, weight) in enumerate(filter.w))) < big"1e-65"
        end
        for solver in (ForwardBackward(), RungeKutta3(), ExactSolver())
            phase = big"0.7"
            response = barotropic_response(phase; solver, filter)
            @test response.W[2, :] - [0, 1] ≈ -im * phase * response.W̃
            @test response.ūₚ[1] ≈ response.W[1, 2]
            @test response.ūₚ⁺[1] ≈ response.W̃[2]
        end
    end
    filter = averaging_filter(NoFilter(); float_type=BigFloat)
    phase = big"0.7"
    response = barotropic_response(phase; solver=ExactSolver(), filter, sampling=ContinuousSampling())
    @test response.W ≈ exp(ComplexF64[0 -0.7im; -0.7im 0]) rtol=1e-14
    @test response.ūₚ[2] ≈ -im * (1 - cos(phase)) / phase
    @test response.ūₚ[3] ≈ -2im * (phase - sin(phase)) / phase^2
    zero_response = barotropic_response(big"0.0"; filter)
    @test zero_response.W == I
    @test zero_response.W̃ == [1, 0]
end

@testset "Stage matrices versus direct updates" begin
    for scheme in (WickerSkamarock3(), ShuOsher3(), ModifiedRungeKutta4(), ModifiedRungeKutta4(a=1 // 4)), solver in (ForwardBackward(), RungeKutta3(), ExactSolver()), count in (1, 4)
        modes = vertical_modes(0.05, count)
        filter = averaging_filter(F3Filter())
        initial = ComplexF64[sin(index) + im * cos(index) for index in 1:2(count+1)]
        forcings = scheme == ShuOsher3() ? (FrozenForcing(), SSPForcing()) : (FrozenForcing(), FinalStageQuadraticForcing(), ProgressiveForcing())
        for forcing in forcings
            matrix = evolution_matrix(0.7, modes; scheme, forcing, solver, filter)
            @test matrix * initial ≈ direct_step(0.7, modes, initial; scheme, forcing, solver, filter) atol=2e-14
        end
    end
    modes = vertical_modes(0.05, 4)
    initial = ComplexF64[sin(index) + im * cos(index) for index in 1:10]
    matrix = evolution_matrix(0.7, modes; forcing=FrozenForcing(), solver=ImplicitSolver())
    @test matrix * initial ≈ direct_step(0.7, modes, initial; forcing=FrozenForcing(), solver=ImplicitSolver())
    filter = averaging_filter(F1Filter())
    history = [initial; reverse(initial)]
    matrix = evolution_matrix(0.7, modes; scheme=AdamsBashforth2(), forcing=FrozenForcing(), solver=ForwardBackward(), filter)
    @test matrix * history ≈ direct_adams_bashforth(0.7, modes, history, filter)
    @test evolution_matrix(0.0, modes) ≈ I atol=1e-13
end

@testset "Eigenvalue perturbation and accuracy" begin
    setprecision(256) do
        modes = vertical_modes(big"0.0001", 4)
        filter = averaging_filter(NoFilter(); float_type=BigFloat)
        matrix = evolution_matrix(big"2.0", modes; forcing=FrozenForcing(), solver=ExactSolver(), filter)
        reference = uncoupled_matrix(big"2.0", modes; solver=ExactSolver(), filter)
        corrections = coupling_corrections(matrix, reference, modes)
        @test corrections.left * corrections.right ≈ I
        @test norm(reference * corrections.right - corrections.right * Diagonal(corrections.eigenvalues)) < big"1e-65"
        @test abs(corrections.coupling - corrections.direct_magnitude - corrections.indirect_magnitude) < abs(corrections.coupling) / 10000
        prediction = frozen_predictions(big"2.0", modes; solver=ExactSolver(), filter)
        @test isapprox(corrections.direct_magnitude, prediction.direct; rtol=big"0.001")
        @test isapprox(corrections.indirect_magnitude, prediction.indirect; rtol=big"0.001")
        @test corrections.overlap > big"0.99"
        modes = vertical_modes(big"0.05", 2)
        for name in ("AB2-SE", "WRK3-SE-F1", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM")
            configuration = scheme_configuration(name)
            filter = averaging_filter(configuration.averaging; float_type=BigFloat)
            errors = [global_error(phase, modes; configuration.scheme, configuration.forcing, configuration.solver, filter, horizon=4).error for phase in (big"0.01", big"0.005", big"0.0025")]
            order = log2(errors[2] / errors[3])
            expected = name in ("AB2-SE", "WRK3-SE-F1", "WRK3-IM") ? 1 : 2
            @test abs(order - expected) < 0.15
        end
    end
end

@testset "Analytical transfers and limits" begin
    setprecision(256) do
        filter = averaging_filter(NoFilter(); float_type=BigFloat)
        for scheme in (WickerSkamarock3(), ShuOsher3())
            errors = BigFloat[]
            for stratification in (big"0.0001", big"0.000025")
                modes = vertical_modes(stratification, 4)
                forcing = FrozenForcing()
                matrix = evolution_matrix(big"0.7", modes; scheme, forcing, solver=ExactSolver(), filter)
                reference = uncoupled_matrix(big"0.7", modes; scheme, solver=ExactSolver(), filter)
                result = coupling_corrections(matrix, reference, modes)
                prediction = frozen_predictions(big"0.7", modes; scheme, solver=ExactSolver(), filter)
                measured = [result.transfers[result.selected, 1:2]; result.transfers[1:2, result.selected]]
                predicted = [prediction.forward; prediction.reverse]
                push!(errors, norm(measured - predicted) / norm(predicted))
            end
            @test 1.8 < errors[1] / errors[2] < 2.2
            modes = vertical_modes(big"0.0001", 4)
            prediction = frozen_predictions(big"0.7", modes; scheme, solver=ImplicitSolver())
            phase = big"0.7"
            coefficient = scheme == WickerSkamarock3() ? (phase^2 + 12) / (12 * (phase^2 + 4)) : (7phase^2 + 6) / (24 * (phase^2 + 1))
            @test prediction.total ≈ -modes.χ[2]^2 * phase^2 * coefficient
        end
        modes = vertical_modes(big"0.0001", 2)
        for scheme in (WickerSkamarock3(), ShuOsher3(), ModifiedRungeKutta4()), solver in (ExactSolver(), RungeKutta3()), phase in (big"0.1", big"0.7", big"2")
            filter = averaging_filter(solver isa ExactSolver ? NoFilter() : F3Filter(); float_type=BigFloat)
            prediction = frozen_predictions(phase, modes; scheme, solver, filter)
            W = barotropic_response(phase; solver, filter).W
            λᵦ, Rᵦ = eigen(W)
            δλ² = sum(prediction.forward .* prediction.reverse ./ (1 .- λᵦ))
            @test prediction.indirect ≈ real(δλ²) rtol=big"1e-50"
            @test prediction.total ≈ prediction.direct + prediction.indirect rtol=big"1e-50"
        end
        filter = averaging_filter(NoFilter(); float_type=BigFloat)
        modes = vertical_modes(big"0.0001", 2)
        for forcing in (FrozenForcing(), FinalStageQuadraticForcing(sampling=ContinuousSampling()), ProgressiveForcing(sampling=ContinuousSampling()))
            phase = big"0.01"
            matrix = evolution_matrix(phase, modes; scheme=ModifiedRungeKutta4(), forcing, solver=ExactSolver(), filter)
            reference = uncoupled_matrix(phase, modes; scheme=ModifiedRungeKutta4(), solver=ExactSolver(), filter)
            result = coupling_corrections(matrix, reference, modes)
            coefficient = result.coupling / (phase * modes.α[2])^4
            expected = forcing == FrozenForcing() ? -1 // 48 : forcing isa FinalStageQuadraticForcing ? 1 // 48 : 0
            @test abs(coefficient - expected) < big"1e-5"
        end
    end
end

@testset "Outer amplification and filter weights" begin
    for ξ in (0, 1 // 10), phase in (big"0", big"0.1", big"1", big"3")
        factor = amplification_factor(phase, AdamsBashforth2(); ξ)
        @test abs(factor^2 - (1 - im * phase * (3 // 2 + ξ)) * factor - im * phase * (1 // 2 + ξ)) < big"1e-65"
    end
    @test amplification_factor(big"0", AdamsBashforth2()) == 1
    @test abs(amplification_factor(big"0.01", AdamsBashforth2()) - 1) < big"0.011"
    for phase in (big"0.1", big"1", big"3")
        @test amplification_factor(phase, WickerSkamarock3()) == amplification_factor(phase, ShuOsher3())
        @test amplification_factor(phase, RungeKutta2()) == 1 - im * phase - phase^2 / 2
        @test amplification_factor(phase, ModifiedRungeKutta4(a=1 // 4)) ≈ sum((-im * phase)^order / factorial(order) for order in 0:4)
    end
    for averaging in (F1Filter(), F2Filter(), F3Filter())
        filter = averaging_filter(averaging, 64; float_type=BigFloat)
        transport = filter.Δτ .* reverse(cumsum(reverse(filter.w)))
        @test sum(transport) ≈ 1
        @test transport - [transport[2:end]; 0] ≈ filter.Δτ .* filter.w
        for phase in (big"0.1", big"2", big"16")
            response = barotropic_response(phase; solver=ExactSolver(), filter)
            factor = sum(weight * exp(-im * phase * index * filter.Δτ) for (index, weight) in enumerate(filter.w))
            @test maximum(abs, eigvals(response.W)) ≈ abs(factor)
        end
    end
end
@testset "Explicit progressive transfers" begin
    setprecision(256) do
        for scheme in (WickerSkamarock3(), ModifiedRungeKutta4(a=23//100), ModifiedRungeKutta4(a=1//4)),
            μ₀ in BigFloat.((2, 6)), propagation in (-1, 1),
            (solver, averaging, sampling) in ((ExactSolver(), NoFilter(), ContinuousSampling()),
                                             (RungeKutta3(), F3Filter(), MidpointSampling()))
            filter = averaging_filter(averaging; float_type=BigFloat)
            errors = map((big"0.0001", big"0.00001")) do ε
                modes = vertical_modes(ε, 4)
                M = evolution_matrix(μ₀, modes; scheme, solver, filter, forcing=ProgressiveForcing(; sampling))
                Mᵤ = uncoupled_matrix(μ₀, modes; scheme, solver, filter)
                result = coupling_corrections(M, Mᵤ, modes; propagation)
                prediction = progressive_predictions(μ₀, modes; scheme, solver, filter, sampling, propagation)
                return max(abs(prediction.direct - result.direct_magnitude),
                           abs(prediction.indirect - result.indirect_magnitude)) / ε^3
            end
            @test errors[2] < errors[1] / 5
        end
        modes = vertical_modes(big"0.01", 4)
        wrk3 = progressive_predictions(big"2.0", modes; scheme=WickerSkamarock3())
        mrk4 = progressive_predictions(big"2.0", modes; scheme=ModifiedRungeKutta4(a=0))
        @test wrk3.direct ≈ mrk4.direct
        @test wrk3.indirect ≈ mrk4.indirect
    end
end

include("rotation.jl")
