using Oceananigans.BuoyancyFormulations: LinearSeawaterBuoyancy

@kernel function _cache_advective_fluxes!(Fⁿ, Fⁿ⁻¹, grid::AbstractGrid, advection, U, buoyancy, C)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        # Save previous advective fluxes
        Fⁿ⁻¹.x[i, j, k] = Fⁿ.x[i, j, k]
        Fⁿ⁻¹.y[i, j, k] = Fⁿ.y[i, j, k]
        Fⁿ⁻¹.z[i, j, k] = Fⁿ.z[i, j, k]

        # Calculate new advective fluxes
        FTx = _advective_tracer_flux_x(i, j, k, grid, advection, U.u, C.T) * σⁿ(i, j, k, grid, f, c, c)
        FTy = _advective_tracer_flux_y(i, j, k, grid, advection, U.v, C.T) * σⁿ(i, j, k, grid, c, f, c)
        FTz = _advective_tracer_flux_z(i, j, k, grid, advection, U.w, C.T) * σⁿ(i, j, k, grid, c, c, f)

        FSx = _advective_tracer_flux_x(i, j, k, grid, advection, U.u, C.S) * σⁿ(i, j, k, grid, f, c, c)
        FSy = _advective_tracer_flux_y(i, j, k, grid, advection, U.v, C.S) * σⁿ(i, j, k, grid, c, f, c)
        FSz = _advective_tracer_flux_z(i, j, k, grid, advection, U.w, C.S) * σⁿ(i, j, k, grid, c, c, f)

        Fⁿ.x[i, j, k] = buoyancy_flux(FTx, FSx, buoyancy)
        Fⁿ.y[i, j, k] = buoyancy_flux(FTy, FSy, buoyancy)
        Fⁿ.z[i, j, k] = buoyancy_flux(FTz, FSz, buoyancy)
    end
end

@kernel function _cache_advective_fluxes!(Fⁿ, grid::AbstractGrid, advection, U, buoyancy, C)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        # Calculate new advective fluxes
        FTx = _advective_tracer_flux_x(i, j, k, grid, advection, U.u, C.T) * σⁿ(i, j, k, grid, f, c, c)
        FTy = _advective_tracer_flux_y(i, j, k, grid, advection, U.v, C.T) * σⁿ(i, j, k, grid, c, f, c)
        FTz = _advective_tracer_flux_z(i, j, k, grid, advection, U.w, C.T) * σⁿ(i, j, k, grid, c, c, f)

        FSx = _advective_tracer_flux_x(i, j, k, grid, advection, U.u, C.S) * σⁿ(i, j, k, grid, f, c, c)
        FSy = _advective_tracer_flux_y(i, j, k, grid, advection, U.v, C.S) * σⁿ(i, j, k, grid, c, f, c)
        FSz = _advective_tracer_flux_z(i, j, k, grid, advection, U.w, C.S) * σⁿ(i, j, k, grid, c, c, f)

        Fⁿ.x[i, j, k] = buoyancy_flux(FTx, FSx, buoyancy)
        Fⁿ.y[i, j, k] = buoyancy_flux(FTy, FSy, buoyancy)
        Fⁿ.z[i, j, k] = buoyancy_flux(FTz, FSz, buoyancy)
    end
end

@inline function buoyancy_flux(FT, FS, buoyancy::LinearSeawaterBuoyancy)
    g = buoyancy.gravitational_acceleration
    α = buoyancy.equation_of_state.thermal_expansion
    β = buoyancy.equation_of_state.haline_contraction
    return g * (α * FT - β * FS)
end
