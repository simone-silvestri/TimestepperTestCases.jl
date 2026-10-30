function substep_response(μₘ, solver)
    W = substep_matrix(μₘ, solver)
    F = W[:, 2] - [0, 1]
    return W, F
end

substep_matrix(μ, ::ForwardBackward) = [one(complex(μ)) -im*μ; -im*μ 1-μ^2]

function substep_matrix(μ, ::ExactSolver)
    J = [0 1; 1 0]
    return cos(μ) * I - im * sin(μ) * J
end

function substep_matrix(μ, ::RungeKutta3)
    J = [0 1; 1 0]
    return I - im * μ * J + (-im * μ * J)^2 / 2 + (-im * μ * J)^3 / 6
end

function substep_matrix(μ, ::ImplicitSolver)
    J = [0 1; 1 0]
    return inv(I + im * μ * J)
end

"""
$(TYPEDSIGNATURES)

Return the filtered barotropic response at barotropic Courant number `μₘ = k√(gH) * stage_duration` in normalized velocity–free-surface coordinates.
`W` is the barotropic evolution matrix W mapping the two initial amplitudes to the final state; `W̃` maps them to the continuity transport.
`ūₚ[p + 1]` and `ūₚ⁺[p + 1]` are the zero-initial-state outputs for unit forcing sᵖ,
where s = t′/stage_duration and p = 0, 1, 2. The superscript ⁺ denotes the manuscript’s dagger continuity transport.
Multiply these outputs by the corresponding forcing coefficients.

`MidpointSampling()` holds the forcing during each substep. `ContinuousSampling()` requires `ExactSolver()` and an unfiltered solve.
Transport is normalized by H√(gH).
"""
function barotropic_response(μₘ; solver=RungeKutta3(), filter=averaging_filter(F3Filter(); float_type=typeof(float(μₘ))), sampling=MidpointSampling())
    T = typeof(float(μₘ))
    if iszero(μₘ)
        return (; W=Matrix{Complex{T}}(I, 2, 2), W̃=Complex{T}[1, 0],
                  ūₚ=zeros(Complex{T}, 3), ūₚ⁺=zeros(Complex{T}, 3))
    end

    W, U = forced_response(μₘ, solver, filter, sampling)

    W̃ = (W[2, :] - T[0, 1]) / (-im * μₘ)
    return (; W, W̃, ūₚ=U[1, :], ūₚ⁺=U[2, :] / (-im * μₘ))
end

function forced_response(μₘ, solver::ExactSolver, filter, ::ContinuousSampling)
    all(iszero, filter.w[1:end-1]) || throw(ArgumentError("Continuous forcing requires an unfiltered solve."))
    T = typeof(float(μₘ))
    W, F₀ = substep_response(μₘ, ExactSolver())
    J = [0 1; 1 0]
    L = -im * μₘ * J
    U = zeros(Complex{T}, 2, 3)
    U[:, 1] = F₀

    for p in 1:2
        if abs(μₘ) <= 1
            term = L[:, 2] / (p + 1)
            series_sum = copy(term)
            for q in 0:precision(T)
                term = L * term / (q + p + 2)
                series_sum += term
                norm(term) < eps(T) * norm(series_sum) && break
            end
            U[:, p + 1] = series_sum
        else
            U[:, p + 1] = p * (L \ U[:, p]) - T[0, 1]
        end
    end
    return W, U
end

function forced_response(μₘ, solver, filter, ::MidpointSampling)
    T     = typeof(float(μₘ))
    Wₛ, F = substep_response(μₘ * filter.Δτ, solver)
    Wˢ    = Matrix{Complex{T}}(I, 2, 2)
    Uˢ    = zeros(Complex{T}, 2, 3)
    W     = zeros(Complex{T}, 2, 2)
    U     = zeros(Complex{T}, 2, 3)

    for (s, w) in enumerate(filter.w)
        τ = (s - T(1 // 2)) * filter.Δτ
        Uˢ = Wₛ * Uˢ + F * transpose(T[1, τ, τ^2])
        Wˢ = Wₛ * Wˢ
        W += w * Wˢ
        U += w * Uˢ
    end

    return W, U
end
