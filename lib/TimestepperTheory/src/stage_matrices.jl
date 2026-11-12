"""
$(TYPEDSIGNATURES)

Return the Section 4 solver, outer scheme, filter, and forcing choices for a manuscript configuration.
"""
function scheme_configuration(name)
    if name == "AB2-SE"
        return (; scheme=AdamsBashforth2(), forcing=FrozenForcing(), solver=ForwardBackward(), averaging=F1Filter())
    elseif name == "WRK3-SE-F1"
        return (; scheme=WickerSkamarock3(), forcing=FrozenForcing(), solver=ForwardBackward(), averaging=F1Filter())
    elseif name == "WRK3-SE" || name == "WRK3-UP"
        return (; scheme=WickerSkamarock3(), forcing=ProgressiveForcing(), solver=RungeKutta3(), averaging=F3Filter())
    elseif name == "SRK3-SE"
        return (; scheme=ShuOsher3(), forcing=FrozenForcing(), solver=RungeKutta3(), averaging=F3Filter())
    elseif name == "MRK4-SE"
        return (; scheme=ModifiedRungeKutta4(), forcing=ProgressiveForcing(), solver=RungeKutta3(), averaging=F3Filter())
    elseif name == "WRK3-IM"
        return (; scheme=WickerSkamarock3(), forcing=FrozenForcing(), solver=ImplicitSolver(), averaging=NoFilter())
    end
    throw(ArgumentError("Unknown manuscript configuration $name."))
end

stage_fractions(::WickerSkamarock3, T)         = T.((1 // 3, 1 // 2, 1))
stage_fractions(scheme::ModifiedRungeKutta4, T) = T.((scheme.a, 1 // 3, 1 // 2, 1))
stage_fractions(::ShuOsher3, T)                = T.((1, 1, 1))

function reconstruction_weights(m, γ)
    T = eltype(γ)
    C = zeros(T, 3, m)
    p = min(m - 1, 2)
    τ = [zero(T); collect(γ[1:m-1])] / γ[m]
    V = [τⱼ^q for τⱼ in τ, q in 0:p]
    C[1:p+1, :] = V \ Matrix{T}(I, m, m)
    return C
end

"""
$(TYPEDSIGNATURES)

Return arrays `(A, B, M)` with A[m] = Aₘ, B[m] = Bₘ, and M[m] = Mₘ.
They satisfy Mₘ = Aₘ + BₘMₘ₋₁ with M₀ = I; earlier-stage reconstructed forcing is included in Aₘ.
`μ₀` is the barotropic Courant number μ₀ = k√(gH)Δt, and `modes` supplies the normalized modal basis in state order [u; h].

`scheme.a` is the first-stage fraction for MRK4. Polynomial forcing uses local time s = t′/stage_duration.
The implicit formulation supports WickerSkamarock3() with FrozenForcing(). SRK3 uses the blended matrices with χᵀχ ≈ 1.
C contains polynomial reconstruction weights, and the rows of 𝓕 multiply yⁿ. K is the frozen-forcing block in Bₘ.
O is the zero block, and M[j + 1] stores the manuscript's Mⱼ.
For the implicit solve, U = [I O] and H = [O I] select the modal amplitudes; ηⁿ and η⁻ multiply yⁿ and yᵐ⁻¹.
"""
function stage_matrices(μ₀, modes; scheme=WickerSkamarock3(), forcing=FinalStageQuadraticForcing(), solver=RungeKutta3(),
                        filter=averaging_filter(F3Filter(); float_type=typeof(float(μ₀))))
    return stage_matrices(solver, scheme, forcing, μ₀, modes; filter)
end

function stage_matrices(solver::Union{ForwardBackward, RungeKutta3, ExactSolver}, scheme, forcing, μ₀, modes; filter)
    T  = typeof(float(μ₀))
    N  = length(modes.α)
    O  = zeros(Complex{T}, N, N)
    M₀ = Matrix{Complex{T}}(I, 2N, 2N)
    
    # Modal weights
    χ,  ψ,  Π  = modes.χ, modes.ψ, modes.Π
    χᵀ, ψᵀ, Πᵀ = transpose(χ), transpose(ψ), transpose(Π)
    𝒫 = I - χ * χᵀ
    
    D = Diagonal(modes.α.^2)
    γ = stage_fractions(scheme, T)
    M = [M₀]
    A = Matrix{Complex{T}}[]
    B = Matrix{Complex{T}}[]

    for m in eachindex(γ)
        μₘ = γ[m] * μ₀
        W, W̃, ūₚ, ūₚ⁺ = barotropic_response(μₘ; solver, filter, sampling=forcing_sampling(forcing))
        
        # Stage A matrix
        Aₘ = [𝒫 + W[1, 1] * χ * χᵀ             W[1, 2] * χ * ψᵀ;
              -im * μₘ * W̃[1] * D * χ * χᵀ     I - im * μₘ * W̃[2] * D * χ * ψᵀ]

        # Stage B matrix
        Bₘ = [O                    -im * μₘ * 𝒫 + W[1, 2] * χ * Πᵀ;
              -im * μₘ * D * 𝒫     -im * μₘ * W̃[2] * D * χ * Πᵀ]

        # Possibly add the forcing reconstruction or blend the stage
        Aₘ, Bₘ = reconstruct_forcing(forcing, m, γ, M; Aₘ, Bₘ, μₘ, χ, D, 𝒫, Πᵀ, W, W̃, ūₚ, ūₚ⁺)
        Aₘ, Bₘ = blend_stage(scheme, m, Aₘ, Bₘ)

        push!(A, Aₘ)
        push!(B, Bₘ)

        # The final stage matric
        Mₘ = Aₘ + Bₘ * M[end]
        push!(M, Mₘ)
    end
    return (; A, B, M=M[2:end])
end

reconstruct_forcing(::FrozenForcing, m, γ, M; Aₘ, Bₘ, kwargs...) = Aₘ, Bₘ

function reconstruct_forcing(::SSPForcing, m, γ, M; Aₘ, Bₘ, μₘ, χ, D, Πᵀ, W, W̃, kwargs...)
    
    # Additional forcing matrix only for the last stage
    if m == 3
        β = (1 // 6, 1 // 6, 2 // 3)
        O = zeros(eltype(Aₘ), length(χ), length(χ))

        K = [O  W[1, 2] * χ * Πᵀ;
             O  -im * μₘ * W̃[2] * D * χ * Πᵀ]
        Aₘ += K * (β[1] * I + β[2] * M[2])
        Bₘ += (β[3] - 1) * K
    end

    return Aₘ, Bₘ
end

function reconstruct_forcing(::FinalStageQuadraticForcing, m, γ, M; Aₘ, Bₘ, kwargs...)

    m < length(γ) && return Aₘ, Bₘ # Reconstruct only the last stage

    return reconstruct_forcing(ProgressiveForcing(), m, γ, M; Aₘ, Bₘ, kwargs...)
end

function reconstruct_forcing(::ProgressiveForcing, m, γ, M; Aₘ, Bₘ, μₘ, χ, D, 𝒫, Πᵀ, ūₚ, ūₚ⁺, kwargs...)
    
    m == 1 && return Aₘ, Bₘ # The first stage is frozen forcing

    N = length(χ)
    O = zeros(eltype(Aₘ), N, N)
    C = reconstruction_weights(m, γ)
    𝓕 = zeros(eltype(Aₘ), 3, 2N)
    for j in 1:m
        𝓕 += C[:, j] * [zeros(eltype(χ), 1, N) Πᵀ] * M[j]
    end
    Aₘ += [χ * transpose(ūₚ) * 𝓕;
           -im * μₘ * D * χ * transpose(ūₚ⁺) * 𝓕]
    Bₘ = [O                    -im * μₘ * 𝒫;
          -im * μₘ * D * 𝒫     O]
    return Aₘ, Bₘ
end

blend_stage(::Union{WickerSkamarock3, ModifiedRungeKutta4}, m, Aₘ, Bₘ) = Aₘ, Bₘ

function blend_stage(::ShuOsher3, m, Aₘ, Bₘ)
    bₘ = (1, 1 // 4, 2 // 3)[m]
    Aₘ = (1 - bₘ) * I + bₘ * (Aₘ - I)
    Bₘ = bₘ * (I + Bₘ)
    return Aₘ, Bₘ
end

function stage_matrices(::ImplicitSolver, scheme::WickerSkamarock3, ::FrozenForcing, μ₀, modes;
                        filter)
    T = typeof(float(μ₀))
    N = length(modes.α)
    O = zeros(Complex{T}, N, N)
    U = [Matrix{Complex{T}}(I, N, N) O]
    H = [O Matrix{Complex{T}}(I, N, N)]
    χ, ψ = modes.χ, modes.ψ
    χᵀ, ψᵀ = transpose(χ), transpose(ψ)
    𝒫 = I - χ * χᵀ
    D = Diagonal(modes.α.^2)
    A = Matrix{Complex{T}}[]
    B = Matrix{Complex{T}}[]
    M = [Matrix{Complex{T}}(I, 2N, 2N)]

    for γₘ in stage_fractions(scheme, T)
        μₘ = γₘ * μ₀
        var"u★" = -im * μₘ * (I - χ * ψᵀ) * H
        ηⁿ = (ψᵀ * H - im * μₘ * χᵀ * U) / (1 + μₘ^2)
        η⁻ = -im * μₘ * χᵀ * var"u★" / (1 + μₘ^2)
        Aₘ = [U - im * μₘ * χ * ηⁿ;
              H - im * μₘ * D * χ * χᵀ * (U - im * μₘ * χ * ηⁿ)]
        Bₘ = [var"u★" - im * μₘ * χ * η⁻;
              -im * μₘ * D * (𝒫 * U + χ * χᵀ * (var"u★" - im * μₘ * χ * η⁻))]
        push!(A, Aₘ)
        push!(B, Bₘ)
        Mₘ = Aₘ + Bₘ * M[end]
        push!(M, Mₘ)
    end
    return (; A, B, M=M[2:end])
end

"""
$(TYPEDSIGNATURES)

Return the complete timestep matrix; AB2 acts on [yⁿ; yⁿ⁻¹] and uses the quasi-AB2 parameter ξ = 1/10 by default.
"""
function evolution_matrix(μ₀, modes; scheme=WickerSkamarock3(), forcing=FinalStageQuadraticForcing(), solver=RungeKutta3(),
                          filter=averaging_filter(F3Filter(); float_type=typeof(float(μ₀))),
                          ξ=1 // 10)
    return evolution_matrix(scheme, μ₀, modes; forcing, solver, filter, ξ)
end

function evolution_matrix(scheme::AbstractRungeKuttaScheme, μ₀, modes; ξ=1 // 10, kwargs...)
    stages = stage_matrices(μ₀, modes; scheme, kwargs...)
    return stages.M[end]
end

function evolution_matrix(::AdamsBashforth2, μ₀, modes; forcing, solver, filter, ξ)
    T = typeof(float(μ₀))
    N = length(modes.α)
    W, W̃, ūₚ, ūₚ⁺ = barotropic_response(μ₀; solver, filter)
    χ = modes.χ
    χᵀ = transpose(χ)
    ψᵀ = transpose(modes.ψ)
    Πᵀ = transpose(modes.Π)
    𝒫 = I - χ * χᵀ
    D = Diagonal(modes.α.^2)

    Aᵤ = [𝒫 + W[1, 1] * χ * χᵀ  W[1, 2] * χ * ψᵀ]
    Aₕ = [W̃[1] * (-im * μ₀ * D * χ) * χᵀ  I + W̃[2] * (-im * μ₀ * D * χ) * ψᵀ]
    A = [Aᵤ; Aₕ]

    Bᵤ = [zeros(Complex{T}, N, N)  -im * μ₀ * 𝒫 + ūₚ[1] * χ * Πᵀ]
    Bₕ = [(-im * μ₀ * D) * 𝒫  ūₚ⁺[1] * (-im * μ₀ * D * χ) * Πᵀ]
    B = [Bᵤ; Bₕ]

    ξ = T(ξ)
    Mⁿ = A + (3 // 2 + ξ) * B
    M⁻ = -(1 // 2 + ξ) * B
    M₀ = Matrix{Complex{T}}(I, 2N, 2N)
    return [Mⁿ M⁻; M₀ zero(M₀)]
end
