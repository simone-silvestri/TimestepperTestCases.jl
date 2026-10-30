function polynomial_filter(Nₛ, r, T)
    τ = collect(1:Nₛ) .* (T(2) / Nₛ)
    τ₀ = T(4 * 8) / (3 * 7)
    w = @. (τ / τ₀)^2 * (1 - (τ / τ₀)^4) - r * τ / τ₀
    resize!(w, findlast(>(0), w))
    w ./= sum(w)
    Δτ = inv(sum(j * wⱼ for (j, wⱼ) in enumerate(w)))
    return (; Δτ, w)
end

"""
$(TYPEDSIGNATURES)

Construct averaging weights for `NoFilter()`, `F1Filter()`, `F2Filter()`, or `F3Filter()` in `float_type`.
Return `(Δτ, w)`, with sample times `j * Δτ` measured in units of the stage duration.
The weights sum to one and their weighted mean sample time is one.

`Nₛ` must be at least eight. For F1–F3 it sets the nominal resolution over [0, 2]; truncation and moment
normalization determine the returned support and spacing. `NoFilter()` uses `Nₛ` samples over [0, 1].
"""
function averaging_filter(filter::AbstractAveragingFilter, Nₛ=48; float_type=Float64)
    Nₛ >= 8 || throw(ArgumentError("At least eight nominal Nₛ are required."))
    return averaging_filter(float_type, filter, Nₛ)
end

function averaging_filter(float_type, ::NoFilter, Nₛ)
    w = zeros(float_type, Nₛ)
    w[end] = 1
    return (; Δτ=inv(float_type(Nₛ)), w)
end

averaging_filter(float_type, ::F1Filter, Nₛ) = polynomial_filter(Nₛ, float_type(18927 // 100000), float_type)

function averaging_filter(float_type, ::F2Filter, Nₛ)
    r₋, r₊ = float_type(18927 // 100000), float_type(3 // 10)
    for n in 1:precision(float_type)+8
        r = (r₋ + r₊) / 2
        filter = polynomial_filter(Nₛ, r, float_type)
        M₂ = sum(wⱼ * (j * filter.Δτ - 1)^2 for (j, wⱼ) in enumerate(filter.w))
        if M₂ > 0
            r₋ = r
        else
            r₊ = r
        end
    end
    return polynomial_filter(Nₛ, (r₋ + r₊) / 2, float_type)
end

function averaging_filter(float_type, ::F3Filter, Nₛ)
    ν  = float_type(595115 // 1000000)
    Δτ = float_type(2) / Nₛ
    τ  = collect(1:Nₛ) .* Δτ

    u = (τ .- float_type(909165 // 1000000)) ./ float_type(864972 // 1000000)
    B = [abs(uᵢ) <= 1 ? (1 - uᵢ^2)^ν * uᵢ^j : zero(uᵢ) for uᵢ in u, j in 0:4]
    M = [sum(B[:, j] .* (τ .- 1).^p) for p in 0:3, j in 1:5] # filter moments
    b = M \ float_type[1, 0, 0, 0]
    w = B * b
    
    resize!(w, findlast(!iszero, w))
    w ./= sum(w)
    Δτ = inv(sum(j * wⱼ for (j, wⱼ) in enumerate(w)))
    return (; Δτ, w)
end
