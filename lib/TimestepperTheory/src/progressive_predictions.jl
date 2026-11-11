"""
$(TYPEDSIGNATURES)

Evaluate the explicit progressive-forcing transfers Vⱼᵦ, Vᵦⱼ, and Vⱼⱼ for WRK3 or MRK4.
Return the direct and barotropic indirect magnitude corrections, retaining terms through order ε³ at fixed μ₀.
The modal coefficients, barotropic solver, averaging filter, and forcing sampling determine the prediction.
"""
function progressive_predictions(μ₀, modes; 
                                 scheme=WickerSkamarock3(), 
                                 selected_mode=1, 
                                 propagation=1,
                                 solver=RungeKutta3(), 
                                 sampling=MidpointSampling(),
                                 filter=averaging_filter(F3Filter(); float_type=typeof(float(μ₀))))
    T = typeof(float(μ₀))
    q = selected_mode + 1
    α, χ, ψ, Π = modes.α[q], modes.χ[q], modes.ψ[q], modes.Π[q]
    χ₀, ψ₀, Π₀ = modes.χ[1], modes.ψ[1], modes.Π[1]
    σ = propagation
    a = scheme isa ModifiedRungeKutta4 ? T(scheme.a) : zero(T)
    γ = [γₘ for γₘ in stage_fractions(scheme, T) if !iszero(γₘ)]
    M = length(γ)
    μ = μ₀ .* γ
    outputs = [barotropic_response(μₘ; solver, filter, sampling) for μₘ in μ]
    W   = getproperty.(outputs, :W)
    W̃   = getproperty.(outputs, :W̃)
    ūₚ  = getproperty.(outputs, :ūₚ)
    ūₚ⁺ = getproperty.(outputs, :ūₚ⁺)
    C = [reconstruction_weights(m, γ) for m in 1:M]
    u = inv(sqrt(T(2)))
    h = σ * α * u
    ū = χ * u

    Πh = zeros(Complex{T}, M)
    Πh[1] = Π * h
    Πh[2] = Π * (h - im * μ[1] * α^2 * u) + Π₀ * (W[1][2, 1] + im * μ[1]) * ū
    for m in 2:M-1
        hᵐ = h - im * μ[m] * α^2 * u - μ[m] * μ[m-1] * α^2 * h
        m == 3 && (hᵐ += im * μ[3] * μ[2] * μ[1] * α^4 * u)
        Πh[m+1] = Π * hᵐ + Π₀ * (W[m][2, 1] + im * μ[m] + μ₀ * μ[m] * γ[m] * ūₚ⁺[m][2]) * ū
    end
    ūᵐ, ū⁺ = progressive_outputs(Πh, C, W, W̃, ūₚ, ūₚ⁺, ū, ψ * h)

    δu = χ * (ūᵐ[M] - ū + μ₀^2 / (2ψ₀) * (ū⁺[M-1] - ū) + im * μ₀ * χ * h * (1 + μ₀^2 / (6ψ₀)))
    δu += χ * α^2 * ū * (a * μ₀^4 / (6ψ₀) + μ₀^2 / 2)
    δu -= μ₀^2 / 2 * α^2 * χ * (ū⁺[M-1] - ū + a * μ₀^2 / 3 * (ū⁺[1] - ū))
    δh = -im * μ₀ * α^2 * χ * (ū⁺[M] - ū + μ₀^2 / 6 * (ū⁺[M-2] - ū))
    δh += μ₀^2 / 2 * α^2 * χ^2 * h * (1 + a * μ₀^2 / 3)
    Vⱼⱼ = (δu + σ * δh / α) / sqrt(T(2))

    δu₀ = (ūᵐ[M] - ū + im * μ₀ * χ * h + μ₀^2 / 2 * α^2 * ū) / χ₀
    δh₀ = (-im * μ₀ * (ū⁺[M] - ū) + μ₀^2 / 2 * χ * h - im * μ₀^3 / 6 * α^2 * ū) / ψ₀
    λᵦ, R = eigen(W[M])
    L = inv(R)
    Vⱼᵦ = zeros(Complex{T}, 2)
    Vᵦⱼ = similar(Vⱼᵦ)
    Wₐ = barotropic_response(a * μ₀; solver, filter, sampling).W
    W⅓ = W[M-2]
    λⱼ = stability_polynomial(-im * σ * μ₀ * α; scheme)

    for b in 1:2
        ūᵦ = χ₀ * R[1, b]
        ηᵦ = ψ₀ * R[2, b]
        Πh[1] = Π₀ * R[2, b]
        for m in 1:M-1
            Πh[m+1] = Π₀ / ψ₀ * (W[m][2, 1] * ūᵦ + W[m][2, 2] * ηᵦ)
        end
        ūᵐ, ū⁺ = progressive_outputs(Πh, C, W, W̃, ūₚ, ūₚ⁺, ūᵦ, ηᵦ)
        ηᵃ = Wₐ[2, 1] * ūᵦ + Wₐ[2, 2] * ηᵦ
        η⅓ = W⅓[2, 1] * ūᵦ + W⅓[2, 2] * ηᵦ
        δu = χ * (ūᵐ[M] - ūᵦ + im * μ₀ * χ₀ / ψ₀ * ηᵦ + μ₀^2 / 2 * χ₀ / ψ₀ * ū⁺[M-1])
        δu -= χ * μ₀^2 / 2 * α^2 * (ū⁺[M-1] - ūᵦ + im * μ₀ / 3 * ηᵃ)
        δh = -im * μ₀ * α^2 * χ * (ū⁺[M] - ūᵦ + im * μ₀ / 2 * η⅓)
        Vⱼᵦ[b] = (δu + σ * δh / α) / sqrt(T(2))
        Vᵦⱼ[b] = L[b, 1] * δu₀ + L[b, 2] * δh₀
    end

    δλ¹ = Vⱼⱼ
    δλ² = sum(Vⱼᵦ .* Vᵦⱼ ./ (λⱼ .- λᵦ))
    radial(δλ) = real(conj(λⱼ) * δλ) / abs(λⱼ)
    direct   = radial(δλ¹)
    indirect = radial(δλ²)
    return (; direct, indirect, total=direct + indirect, forward=Vⱼᵦ, reverse=Vᵦⱼ, δλ¹, δλ²)
end

function progressive_outputs(Πh, C, W, W̃, ūₚ, ūₚ⁺, ū, η)
    ūᵐ = similar(Πh)
    ū⁺ = similar(Πh)
    for m in eachindex(W)
        F = C[m] * Πh[1:m]
        ūᵐ[m] = W[m][1, 1] * ū + W[m][1, 2] * η + sum(ūₚ[m] .* F)
        ū⁺[m] = W̃[m][1]    * ū + W̃[m][2]    * η + sum(ūₚ⁺[m] .* F)
    end
    return ūᵐ, ū⁺
end
