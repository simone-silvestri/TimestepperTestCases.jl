"""
$(TYPEDSIGNATURES)

Return the uncoupled evolution matrix: the barotropic evolution matrix W for mode zero and the outer stability polynomial for each baroclinic mode.
"""
function uncoupled_matrix(μ₀, modes; scheme=WickerSkamarock3(), solver=RungeKutta3(),
                          filter=averaging_filter(F3Filter(); float_type=typeof(float(μ₀))), sampling=MidpointSampling(), a=23 // 100)
    N = length(modes.α)
    M = zeros(Complex{typeof(float(μ₀))}, 2N, 2N)
    W = barotropic_window(μ₀, solver; filter, sampling)
    M[[1, N + 1], [1, N + 1]] = W
    for index in 2:N
        Lₓ = -im * μ₀ * [0 1; modes.α[index]^2 0]
        M[[index, N + index], [index, N + index]] = stability_polynomial(Lₓ; scheme, a)
    end
    return M
end

"""
$(TYPEDSIGNATURES)

Return the coupling coefficients Vⱼₗ, direct and indirect eigenvalue corrections, magnitude corrections, and full-matrix growth.
`selected_mode=1` selects the first baroclinic mode; `propagation` is +1 or −1 and selects the sign in (1, ±α).
`direct` is the first-order eigenvalue shift; `indirect` sums second-order corrections through the two barotropic
eigenvalues, and `other` sums the remaining indirect corrections. Their `_magnitude` fields give
Re(conj(λⱼ) δλ)/|λⱼ|, the linearized change in the magnitude of the amplification factor. `growth` is |λ| − 1 and `coupling` is |λ| − |λ₀| for the full eigenmode selected by
maximum energy-scaled eigenvector overlap. `overlap` reports that match; exact eigenvalue degeneracies are singular.
"""
function coupling_corrections(M, Mᵤ, modes; selected_mode=1, propagation=1)
    N = length(modes.α)
    T = eltype(modes.α)
    R = zeros(Complex{T}, 2N, 2N)
    L = zeros(Complex{T}, 2N, 2N)
    λ = zeros(Complex{T}, 2N)

    external_indices = [1, N + 1]
    external = eigen(Mᵤ[external_indices, external_indices])
    R[external_indices, 1:2] = external.vectors
    L[1:2, external_indices] = inv(external.vectors)
    λ[1:2] = external.values

    for q in 1:N-1, (d, σ) in enumerate((1, -1))
        k = 2q + d
        indices = [q + 1, N + q + 1]
        α = modes.α[q + 1]
        R[indices, k] = T[1, σ * α] / sqrt(T(2))
        L[k, indices] = T[1, σ / α] / sqrt(T(2))
        λ[k] = (L[k:k, indices] * Mᵤ[indices, indices] * R[indices, k])[1]
    end

    V = L * (M - Mᵤ) * R
    j = 2selected_mode + (propagation == 1 ? 1 : 2)
    δλ¹ = V[j, j]
    δλ² = zero(δλ¹)
    other = zero(δλ¹)
    for l in eachindex(λ)
        l == j && continue
        λ[j] == λ[l] && throw(DomainError(λ[j], "The perturbation expansion requires distinct eigenvalues."))
        δ = V[j, l] * V[l, j] / (λ[j] - λ[l])
        if l <= 2
            δλ² += δ
        else
            other += δ
        end
    end

    full = eigen(M)
    scaling = Diagonal([ones(T, N); inv.(modes.α)])
    target = scaling * R[:, j]
    overlaps = map(eachcol(full.vectors)) do r
        scaled_vector = scaling * r
        return abs(dot(target, scaled_vector)) / (norm(target) * norm(scaled_vector))
    end
    λ̂ = full.values[argmax(overlaps)]
    λⱼ = λ[j]
    radial(δ) = real(conj(λⱼ) * δ) / abs(λⱼ)
    return (; direct=δλ¹, indirect=δλ², other, transfers=V, eigenvalues=λ, right=R, left=L, selected=j,
              direct_magnitude=radial(δλ¹), indirect_magnitude=radial(δλ²), other_magnitude=radial(other),
              reference_eigenvalue=λⱼ, actual_eigenvalue=λ̂, overlap=maximum(overlaps),
              damping=1 - abs(λⱼ), growth=abs(λ̂) - 1, coupling=abs(λ̂) - abs(λⱼ),
              spectral_growth=maximum(abs, full.values) - 1)
end

"""
$(TYPEDSIGNATURES)

Evaluate the leading weak-stratification frozen-forcing predictions for WRK3, MRK4, or SRK3, including the exact modal weight χ².
"""
function frozen_predictions(μ₀, modes; scheme=WickerSkamarock3(), selected_mode=1, solver=RungeKutta3(),
                            filter=averaging_filter(F3Filter(); float_type=typeof(float(μ₀))))
    W  = barotropic_window(μ₀, solver; filter)
    W½ = barotropic_window(μ₀ / 2, solver; filter)
    δλ¹, total = frozen_corrections(μ₀, scheme, W, W½)

    χ = modes.χ[selected_mode + 1]
    χ² = χ^2
    external = eigen(W)
    L = inv(external.vectors)
    forward = zeros(Complex{typeof(float(μ₀))}, 2)
    reverse = similar(forward)
    for l in 1:2
        λₗ = external.values[l]
        rₗ = external.vectors[:, l]
        ℓₗ = L[l, :]
        reverse[l] = χ / sqrt(typeof(float(μ₀))(2)) *
                           ((λₗ - 1) * ℓₗ[1] + im * μ₀ * (1 + λₗ) * ℓₗ[2] / 2)
        pressure_response = pressure_transfer(scheme, λₗ, rₗ, W½)
        forward[l] = χ / sqrt(typeof(float(μ₀))(2)) * ((λₗ - 1) * rₗ[1] + im * μ₀ * pressure_response)
    end
    return (; direct=χ² * real(δλ¹), indirect=χ² * real(total - δλ¹), total=χ² * real(total), forward, reverse)
end

"""
$(TYPEDSIGNATURES)

Return the fixed-horizon operator error in the Euclidean norm of [u; h/α], with exact initial history for AB2.
`μ₀` and `horizon` are measured in dimensionless time k√(gH)t. The phase is adjusted so an integer
number of steps reaches the horizon exactly. Return `(phase, error)`, including the adjusted phase.
Solver and forcing keywords are passed to `evolution_matrix`; the final operator norm is evaluated in Float64.
"""
function global_error(μ₀, modes; horizon=64, scheme=WickerSkamarock3(), kwargs...)
    steps = max(1, round(Int, horizon / μ₀))
    μ = typeof(float(μ₀))(horizon) / steps
    M = evolution_matrix(μ, modes; scheme, kwargs...)
    Mⁿ = propagated_state(M, steps, scheme, μ, modes)
    scaling = Diagonal([ones(eltype(modes.α), length(modes.α)); inv.(modes.α)])
    E = scaling * (Mⁿ - exact_evolution(typeof(float(μ₀))(horizon), modes)) / scaling
    error = opnorm(ComplexF64.(E))
    return (; phase=μ, error)
end

barotropic_window(μ₀, solver::ImplicitSolver; kwargs...) = substep_matrix(μ₀, solver)

function barotropic_window(μ₀, solver::Union{ForwardBackward, RungeKutta3, ExactSolver}; kwargs...)
    response = barotropic_response(μ₀; solver, kwargs...)
    return response.W
end

function frozen_corrections(μ₀, ::Union{WickerSkamarock3, ModifiedRungeKutta4}, W, W½)
    S = inv(I - W)
    δλ¹ = (W[1, 1] - 1 + im * μ₀ * W[1, 2] / 2 + im * μ₀ * W½[2, 1]
              - μ₀^2 * (W½[2, 2] - 1) / 3 - μ₀^2 / 2) / 2
    total = μ₀^2 * (W½[2, 2] - 1 - 6 * (W½ * S)[2, 2]) / 12
    return δλ¹, total
end

function frozen_corrections(μ₀, ::ShuOsher3, W, W½)
    S = inv(I - W)
    δλ¹ = (W[1, 1] - 1 + im * μ₀ * (W[1, 2] + W[2, 1]) / 2
              - μ₀^2 * (2 + W[2, 2]) / 6) / 2
    total = μ₀^2 * (5 + W[2, 2] - 12 * S[2, 2]) / 24
    return δλ¹, total
end

pressure_transfer(::ShuOsher3, λₗ, r, W½) = (1 + λₗ) * r[2] / 2
pressure_transfer(::Union{WickerSkamarock3, ModifiedRungeKutta4}, λₗ, r, W½) = sum(W½[2, :] .* r)

propagated_state(M, steps, ::AbstractRungeKuttaScheme, μ₀, modes) = M^steps

function propagated_state(M, steps, ::AdamsBashforth2, μ₀, modes)
    dimension = 2length(modes.α)
    history = [Matrix{eltype(M)}(I, dimension, dimension); exact_evolution(-μ₀, modes)]
    return (M^steps * history)[1:dimension, :]
end
