using Oceananigans.BuoyancyFormulations: LinearSeawaterBuoyancy, BuoyancyForce

@inline b★(i, j, k, grid, bⁿ⁺¹, bⁿ) = @inbounds (bⁿ⁺¹[i, j, k] + bⁿ[i, j, k]) / 2
@inline b²(i, j, k, grid, bⁿ⁺¹, bⁿ) = @inbounds bⁿ⁺¹[i, j, k] * bⁿ[i, j, k]

@kernel function _assemble_advective_dissipation!(P, grid, C₁, C₂, Fⁿ, Fⁿ⁻¹, Uⁿ, Uⁿ⁻¹, bⁿ⁺¹, bⁿ)
    i, j, k = @index(Global, NTuple)

    δˣb★ = δxᶠᶜᶜ(i, j, k, grid, b★, bⁿ⁺¹, bⁿ)
    δˣb² = δxᶠᶜᶜ(i, j, k, grid, b², bⁿ⁺¹, bⁿ)

    δʸb★ = δyᶜᶠᶜ(i, j, k, grid, b★, bⁿ⁺¹, bⁿ)
    δʸb² = δyᶜᶠᶜ(i, j, k, grid, b², bⁿ⁺¹, bⁿ)

    δᶻb★ = δzᶜᶜᶠ(i, j, k, grid, b★, bⁿ⁺¹, bⁿ)
    δᶻb² = δzᶜᶜᶠ(i, j, k, grid, b², bⁿ⁺¹, bⁿ)

    @inbounds begin
        fx = C₁ * Fⁿ.x[i, j, k] - C₂ * Fⁿ⁻¹.x[i, j, k]
        fy = C₁ * Fⁿ.y[i, j, k] - C₂ * Fⁿ⁻¹.y[i, j, k]
        fz = C₁ * Fⁿ.z[i, j, k] - C₂ * Fⁿ⁻¹.z[i, j, k]

        u = C₁ * Uⁿ.u[i, j, k] - C₂ * Uⁿ⁻¹.u[i, j, k]
        v = C₁ * Uⁿ.v[i, j, k] - C₂ * Uⁿ⁻¹.v[i, j, k]
        w = C₁ * Uⁿ.w[i, j, k] - C₂ * Uⁿ⁻¹.w[i, j, k]

        P.x[i, j, k] = 2 * δˣb★ * fx - δˣb² * u
        P.y[i, j, k] = 2 * δʸb★ * fy - δʸb² * v
        P.z[i, j, k] = 2 * δᶻb★ * fz - δᶻb² * w
    end
end

#####
##### Fluxes and transports, area-weighted, on the grid the stage computed its tendencies on.
##### The vertical face area does not stretch; the horizontal ones are taken back to the stretching `σᵐ⁻¹`.
#####

@inline horizontal_stretching_ratios(i, j, k, grid, σᵐ⁻¹) =
    @inbounds (σᵐ⁻¹.x[i, j, 1] / σⁿ(i, j, k, grid, f, c, c), σᵐ⁻¹.y[i, j, 1] / σⁿ(i, j, k, grid, c, f, c))

# Adds the contribution `a F` of one tracer flux `F` to the buoyancy flux already in `Fᵇ`, weighted by `keep`
@kernel function _accumulate_advective_fluxes!(Fᵇ, grid::AbstractGrid, advection, U, tracer, σᵐ⁻¹, a, keep)
    i, j, k = @index(Global, NTuple)

    rˣ, rʸ = horizontal_stretching_ratios(i, j, k, grid, σᵐ⁻¹)

    Fx = _advective_tracer_flux_x(i, j, k, grid, advection, U.u, tracer) * rˣ
    Fy = _advective_tracer_flux_y(i, j, k, grid, advection, U.v, tracer) * rʸ
    Fz = _advective_tracer_flux_z(i, j, k, grid, advection, U.w, tracer)

    @inbounds begin
        Fᵇ.x[i, j, k] = keep * Fᵇ.x[i, j, k] + a * Fx
        Fᵇ.y[i, j, k] = keep * Fᵇ.y[i, j, k] + a * Fy
        Fᵇ.z[i, j, k] = keep * Fᵇ.z[i, j, k] + a * Fz
    end
end

@kernel function _accumulate_transport!(U★, grid, U, σᵐ⁻¹, weight, keep)
    i, j, k = @index(Global, NTuple)

    rˣ, rʸ = horizontal_stretching_ratios(i, j, k, grid, σᵐ⁻¹)

    @inbounds begin
        U★.u[i, j, k] = keep * U★.u[i, j, k] + weight * U.u[i, j, k] * Axᶠᶜᶜ(i, j, k, grid) * rˣ
        U★.v[i, j, k] = keep * U★.v[i, j, k] + weight * U.v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid) * rʸ
        U★.w[i, j, k] = keep * U★.w[i, j, k] + weight * U.w[i, j, k] * Azᶜᶜᶠ(i, j, k, grid)
    end
end

@kernel function _store_horizontal_stretching!(σ, grid)
    i, j = @index(Global, NTuple)
    @inbounds begin
        σ.x[i, j, 1] = σⁿ(i, j, 1, grid, f, c, c)
        σ.y[i, j, 1] = σⁿ(i, j, 1, grid, c, f, c)
    end
end

@inline buoyancy_coefficients(bf::BuoyancyForce) = buoyancy_coefficients(bf.formulation)

# The weights of the temperature and salinity fluxes in the buoyancy flux g (α Fᵀ - β Fˢ)
@inline function buoyancy_coefficients(buoyancy::LinearSeawaterBuoyancy)
    g = buoyancy.gravitational_acceleration
    α = buoyancy.equation_of_state.thermal_expansion
    β = buoyancy.equation_of_state.haline_contraction
    return (T = g * α, S = - g * β)
end
