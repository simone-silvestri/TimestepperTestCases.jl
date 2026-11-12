abstract type AbstractRotationCoupling end
struct FrozenRotation <: AbstractRotationCoupling end
struct ProgressivePressure <: AbstractRotationCoupling end
struct ProgressiveRotation <: AbstractRotationCoupling end
struct SSPRotation <: AbstractRotationCoupling end

function rotation_reconstruction(m, c, ::Union{FrozenRotation, SSPRotation})
    C = zeros(3, m)
    C[1, m] = 1
    return C
end

rotation_reconstruction(m, c, ::Union{ProgressivePressure, ProgressiveRotation}) = reconstruction_weights(m, c)

function rotation_outputs(μ, c, filter)
    W, U = forced_response(c * μ, RungeKutta3(), filter, MidpointSampling())
    J  = zeros(3)
    Jˢ = zeros(3)
    for (s, w) in enumerate(filter.w)
        τ = (s - 1/2) * filter.Δτ
        Jˢ += filter.Δτ .* [1, τ, τ^2]
        J += w * Jˢ
    end
    return (; W, U, J=c * J)
end

"""
    rotating_evolution_matrix(μ, ρ, modes, coupling; scheme=WickerSkamarock3(), filter=averaging_filter(F3Filter(), 96), outputs=nothing)

Return the S8 timestep matrix for `(u, v, h)`, with gravity-wave phase `μ = k√(gH)Δt` and inertial phase `ρ = fΔt`.
`FrozenRotation`, `ProgressivePressure`, and `ProgressiveRotation` specify frozen or progressive forcing for the restarted `scheme`.
`SSPRotation` holds pressure and Coriolis at their preceding-stage values.
The filter applies to all barotropic outputs. Prescribed polynomial forcing uses midpoint sampling.
"""
function rotating_evolution_matrix(μ, ρ, modes, coupling; scheme=WickerSkamarock3(), filter=averaging_filter(F3Filter(), 96), outputs=nothing)
    N = length(modes.α)
    χ, ψ, Π = modes.χ, modes.ψ, modes.Π
    D = Diagonal(modes.α.^2)
    P = I - χ * transpose(χ)
    E = Matrix{ComplexF64}(I, 3N, 3N)
    U, V, H = E[1:N, :], E[N+1:2N, :], E[2N+1:3N, :]
    c = coupling isa SSPRotation ? (1., 1., 1.) : stage_fractions(scheme, Float64)
    M = [E]
    for m in eachindex(c)
        W, B, J = isnothing(outputs) ? rotation_outputs(μ, c[m], filter) : outputs[m]
        C  = rotation_reconstruction(m, c, coupling)
        Cᶠ = coupling isa ProgressivePressure ? rotation_reconstruction(m, c, FrozenRotation()) : C
        F  = zeros(ComplexF64, 3, 3N)
        G  = zero(F)
        for j in 1:m
            F += C[:, j] * transpose(Π) * H * M[j]
            F += im * ρ / μ * Cᶠ[:, j] * transpose(χ) * V * M[j]
            G -= ρ * Cᶠ[:, j] * transpose(χ) * U * M[j]
        end
        X = W * [transpose(χ) * U; transpose(ψ) * H] + B * F
        v̄ = transpose(χ) * V + transpose(J) * G
        u = P * (U + c[m] * ρ * V * M[end] - im * c[m] * μ * H * M[end]) + χ * X[1:1, :]
        v = P * (V - c[m] * ρ * U * M[end]) + χ * v̄
        u⁺ = (X[2:2, :] - transpose(ψ) * H) / (-im * c[m] * μ)
        h = H - im * c[m] * μ * D * (P * U * M[end] + χ * u⁺)
        A = [u; v; h]
        if coupling isa SSPRotation
            b = (1., 1/4, 2/3)[m]
            A = (1-b) * E + b * (M[end] + A - E)
        end
        push!(M, A)
    end
    return M[end]
end

