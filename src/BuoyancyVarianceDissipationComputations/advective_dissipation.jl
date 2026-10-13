using Oceananigans.BuoyancyFormulations: LinearSeawaterBuoyancy, BuoyancyForce

# Adds the contribution `a F` of one tracer flux `F` to the buoyancy flux already in `Fᵇ`, weighted by `keep`
@kernel function _accumulate_advective_fluxes!(Fᵇ, grid::AbstractGrid, advection, U, tracer, σ, a, keep)
    i, j, k = @index(Global, NTuple)

    Fx = _advective_tracer_flux_x(i, j, k, grid, advection, U.u, tracer) * σ(i, j, k, grid, f, c, c)
    Fy = _advective_tracer_flux_y(i, j, k, grid, advection, U.v, tracer) * σ(i, j, k, grid, c, f, c)
    Fz = _advective_tracer_flux_z(i, j, k, grid, advection, U.w, tracer) * σ(i, j, k, grid, c, c, f)

    @inbounds begin
        Fᵇ.x[i, j, k] = keep * Fᵇ.x[i, j, k] + a * Fx
        Fᵇ.y[i, j, k] = keep * Fᵇ.y[i, j, k] + a * Fy
        Fᵇ.z[i, j, k] = keep * Fᵇ.z[i, j, k] + a * Fz
    end
end

@kernel function _store_previous_advective_fluxes!(Fⁿ⁻¹, Fⁿ)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        Fⁿ⁻¹.x[i, j, k] = Fⁿ.x[i, j, k]
        Fⁿ⁻¹.y[i, j, k] = Fⁿ.y[i, j, k]
        Fⁿ⁻¹.z[i, j, k] = Fⁿ.z[i, j, k]
    end
end

@inline unit_scaling(i, j, k, grid, ℓx, ℓy, ℓz) = one(eltype(grid))

@inline buoyancy_coefficients(bf::BuoyancyForce) = buoyancy_coefficients(bf.formulation)

# The weights of the temperature and salinity fluxes in the buoyancy flux g (α Fᵀ - β Fˢ)
@inline function buoyancy_coefficients(buoyancy::LinearSeawaterBuoyancy)
    g = buoyancy.gravitational_acceleration
    α = buoyancy.equation_of_state.thermal_expansion
    β = buoyancy.equation_of_state.haline_contraction
    return (T = g * α, S = - g * β)
end
