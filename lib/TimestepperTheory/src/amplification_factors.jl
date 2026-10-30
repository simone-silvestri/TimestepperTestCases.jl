"""
$(TYPEDSIGNATURES)

Return the amplification factor for unit-speed Fourier advection at phase θ = kΔt. Adams–Bashforth uses the physical
root that approaches one as θ approaches zero; `ξ` is its over-extrapolation parameter.
"""
function amplification_factor(θ, scheme::AbstractRungeKuttaScheme; a=23 // 100)
    return stability_polynomial(-im * θ; scheme, a)
end

function amplification_factor(θ, ::AdamsBashforth2; ξ=1 // 10)
    d = 1 - im * (1 - 2ξ) * θ - (3 // 2 + ξ)^2 * θ^2
    return (1 - im * (3 // 2 + ξ) * θ + sqrt(complex(d))) / 2
end

stability_polynomial(z, ::RungeKutta2, a) = one(z) + z + z^2 / 2
