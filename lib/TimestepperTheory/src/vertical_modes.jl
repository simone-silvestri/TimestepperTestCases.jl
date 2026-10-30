"""
$(TYPEDSIGNATURES)

Return the constant-stratification modes for `ε = N²H/g > 0` and baroclinic mode count `Q`.
Index `q + 1` represents mode q; index one is the barotropic mode. The floating type follows `ε`.

The returned fields are `x` (x solving x tan(x) = ε), `A` (normalization), `α` (α = c/√(gH)),
`ψ` (ψ = A cos(x)), `χ` (χ = A sin(x)/x), and `Π` (Π = ψ(α² − 1)).
The modes have unit depth-averaged squared amplitude; finite truncations retain their unmodified averages.
"""
function vertical_modes(ε, Q)
    ε > 0 || throw(ArgumentError("Stratification must be positive."))
    Q >= 0 || throw(ArgumentError("The baroclinic mode count must be nonnegative."))
    T = typeof(float(ε))
    x = zeros(T, Q + 1)
    for q in 0:Q
        a = q * T(π)
        b = a + T(π) / 2
        for iteration in 1:precision(T)+16
            c = (a + b) / 2
            if c * tan(c) < ε
                a = c
            else
                b = c
            end
        end
        x[q + 1] = (a + b) / 2
    end
    A = @. inv(sqrt(1 / 2 + sin(2x) / (4x)))
    α = sqrt(ε) ./ x
    ψ = A .* cos.(x)
    χ = A .* sin.(x) ./ x
    Π = ψ .* (α.^2 .- 1)
    return (; x, A, α, ψ, χ, Π)
end

"""
$(TYPEDSIGNATURES)

Return the exact evolution matrix over barotropic Courant number `μ₀ = k√(gH)Δt` for the supplied `modes`.
The state is [u₀, …, uQ, h₀, …, hQ], with velocity normalized by √(gH) and total pressure by gH.
The returned complex matrix has size 2(Q + 1) × 2(Q + 1).
"""
function exact_evolution(μ₀, modes)
    N = length(modes.α)
    M = zeros(Complex{typeof(float(μ₀))}, 2N, 2N)
    for q in 1:N
        α = modes.α[q]
        θ = μ₀ * α
        M[q, q] = M[N + q, N + q] = cos(θ)
        M[q, N + q] = -im * sin(θ) / α
        M[N + q, q] = -im * sin(θ) * α
    end
    return M
end
