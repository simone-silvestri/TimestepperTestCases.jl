using TimestepperTheory
using TimestepperTheory: forced_response, rotation_outputs, stage_fractions
using LinearAlgebra
using Serialization
using SHA
using Test

BLAS.set_num_threads(1)

function scan_cases()
    cases = [
        (; name="WRK3-SE",                   scheme=WickerSkamarock3(),           coupling=ProgressiveRotation()),
        (; name="SRK3-SE",                   scheme=ShuOsher3(),                  coupling=SSPRotation()),
        (; name="MRK4-SE",                   scheme=ModifiedRungeKutta4(a=0.23),  coupling=ProgressiveRotation()),
        (; name="WRK3-frozen-F3",            scheme=WickerSkamarock3(),           coupling=FrozenRotation()),
        (; name="WRK3-progressive-pressure", scheme=WickerSkamarock3(),           coupling=ProgressivePressure()),
        (; name="MRK4-frozen-F3",            scheme=ModifiedRungeKutta4(a=0.23),  coupling=FrozenRotation()),
        (; name="MRK4-progressive-pressure", scheme=ModifiedRungeKutta4(a=0.23),  coupling=ProgressivePressure()),
    ]

    return [(; case..., solver=RungeKutta3(), averaging=F3Filter(), substeps=96) for case in cases]
end

function stage_outputs(μ, case, filter)
    c = stage_fractions(case.scheme, Float64)
    return map(c) do γ
        W, U = forced_response(γ * μ, case.solver, filter, MidpointSampling())
        J = rotation_outputs(μ, γ, filter).J
        (; W, U, J)
    end
end

function amplification_matrix(μ, ρ, modes, case, filter, outputs=stage_outputs(μ, case, filter))
    return rotating_evolution_matrix(μ, ρ, modes, case.coupling; scheme=case.scheme, filter, outputs)
end

function outer_limit(scheme)
    R(θ) = stability_polynomial(im * θ; scheme)
    a, b = 1.0, 4.0
    for n in 1:55
        c = (a + b) / 2
        abs(R(c)) > 1 ? (b = c) : (a = c)
    end
    return (a + b) / 2
end

function validate_cases(cases)
    modes   = vertical_modes(0.1, 8)
    N       = length(modes.α)
    indices = [collect(1:N); collect((2N+1):3N)]
    @testset "Zero-rotation and sign checks" begin
        for case in cases
            filter  = averaging_filter(case.averaging, case.substeps)
            forcing = case.coupling isa Union{FrozenRotation, SSPRotation} ? FrozenForcing() : ProgressiveForcing()
            A       = amplification_matrix(2.0, 0.0, modes, case, filter)
            B       = evolution_matrix(2.0, modes; scheme=case.scheme, forcing, solver=case.solver, filter)
            @test A[indices, indices] ≈ B atol=5e-12
            @test A[(N+1):2N, (N+1):2N] ≈ I atol=5e-12
            for μ in (0.03, 3.0, 20.0)
                A = amplification_matrix(μ, 0.4, modes, case, filter)
                B = amplification_matrix(μ, -0.4, modes, case, filter)
                @test sort(abs.(eigvals(A))) ≈ sort(abs.(eigvals(B))) atol=5e-11
            end
        end
    end
end

function run_scan(directory; cases=scan_cases())
    mkpath(directory)
    validate_cases(cases)
    ρ       = collect(range(0.0, 2.0; length=201))
    results = NamedTuple[]
    open(joinpath(directory, "source_hashes.txt"), "w") do io
        root = dirname(pathof(TimestepperTheory))
        for file in sort(readdir(root; join=true))
            endswith(file, ".jl") && println(io, bytes2hex(sha256(read(file))), "  ", file)
        end
    end
    for ε in (0.01, 0.1), case in cases
        modes = vertical_modes(ε, 8)
        θ★ = outer_limit(case.scheme)
        μ = sort(unique([10.0 .^ range(-4.0, -1.0; length=13); collect(range(0.1, 6.0; length=120));
                         collect(range(6.0, 1.05θ★/modes.α[2]; length=181))]))
        θ       = modes.α[2] .* μ
        filter  = averaging_filter(case.averaging, case.substeps)
        growth  = zeros(length(μ), length(ρ))
        elapsed = @elapsed Threads.@threads for i in eachindex(μ)
            outputs = stage_outputs(μ[i], case, filter)
            for j in eachindex(ρ)
                growth[i, j] = maximum(abs, eigvals(amplification_matrix(μ[i], ρ[j], modes, case, filter, outputs))) - 1
            end
        end
        result = (; case, ε, modes, θ★, μ, θ, ρ, growth)
        push!(results, result)
        serialize(joinpath(directory, "scan.pending.jls"), results)
        mv(joinpath(directory, "scan.pending.jls"), joinpath(directory, "scan.jls"); force=true)
        open(joinpath(directory, "$(case.name)-epsilon$(ε).csv"), "w") do io
            println(io, "mu0,theta1,fdt,growth")
            for j in eachindex(ρ), i in eachindex(μ)
                println(io, μ[i], ',', θ[i], ',', ρ[j], ',', growth[i, j])
            end
        end
        println("COMPLETED ", case.name, " epsilon=", ε, " substeps=", case.substeps, " seconds=",
                round(elapsed; digits=1))
        flush(stdout)
    end
    return results
end
