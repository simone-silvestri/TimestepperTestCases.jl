"""
$(TYPEDSIGNATURES)

Evaluate the outer Runge–Kutta stability polynomial at a scalar or square matrix argument.
"""
stability_polynomial(z; scheme=WickerSkamarock3())            = stability_polynomial(z, scheme)
stability_polynomial(z, ::Union{WickerSkamarock3, ShuOsher3}) = one(z) + z + z^2 / 2 + z^3 / 6
stability_polynomial(z, scheme::ModifiedRungeKutta4)          = one(z) + z + z^2 / 2 + z^3 / 6 + scheme.a * z^4 / 6

"""
$(TYPEDSIGNATURES)

Return the amplification factor for unit-speed Fourier advection at phase θ = kΔt. Adams–Bashforth uses the physical
root that approaches one as θ approaches zero; `ξ` is its over-extrapolation parameter.
"""
amplification_factor(θ, scheme::AbstractRungeKuttaScheme) = stability_polynomial(-im * θ; scheme)

function amplification_factor(θ, ::AdamsBashforth2; ξ=1 // 10)
    d = 1 - im * (1 - 2ξ) * θ - (3 // 2 + ξ)^2 * θ^2
    return (1 - im * (3 // 2 + ξ) * θ + sqrt(complex(d))) / 2
end

stability_polynomial(z, ::RungeKutta2) = one(z) + z + z^2 / 2
