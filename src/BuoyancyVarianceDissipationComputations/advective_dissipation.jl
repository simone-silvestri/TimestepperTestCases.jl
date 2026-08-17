using Oceananigans.BuoyancyFormulations: LinearSeawaterBuoyancy, BuoyancyForce

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

@inline buoyancy_flux(FT, FS, bf::BuoyancyForce) = buoyancy_flux(FT, FS, bf.formulation)

@inline function buoyancy_flux(FT, FS, buoyancy::LinearSeawaterBuoyancy)
    g = buoyancy.gravitational_acceleration
    α = buoyancy.equation_of_state.thermal_expansion
    β = buoyancy.equation_of_state.haline_contraction
    return g * (α * FT - β * FS)
end

# Strong-stability-preserving accumulation of the buoyancy flux. Contrarily to the kernels above, the fluxes
# are summed RAW, without the grid scaling σⁿ: the Shu-Osher path leaves the σ cache at unity, so applying σ
# here as well would scale the budget twice on a moving grid.
@kernel function _accumulate_ssp_advective_fluxes!(F★, grid::AbstractGrid, advection, U, buoyancy, C, β, keep)
    i, j, k = @index(Global, NTuple)

    FTx = _advective_tracer_flux_x(i, j, k, grid, advection, U.u, C.T)
    FTy = _advective_tracer_flux_y(i, j, k, grid, advection, U.v, C.T)
    FTz = _advective_tracer_flux_z(i, j, k, grid, advection, U.w, C.T)

    FSx = _advective_tracer_flux_x(i, j, k, grid, advection, U.u, C.S)
    FSy = _advective_tracer_flux_y(i, j, k, grid, advection, U.v, C.S)
    FSz = _advective_tracer_flux_z(i, j, k, grid, advection, U.w, C.S)

    @inbounds begin
        F★.x[i, j, k] = keep * F★.x[i, j, k] + β * buoyancy_flux(FTx, FSx, buoyancy)
        F★.y[i, j, k] = keep * F★.y[i, j, k] + β * buoyancy_flux(FTy, FSy, buoyancy)
        F★.z[i, j, k] = keep * F★.z[i, j, k] + β * buoyancy_flux(FTz, FSz, buoyancy)
    end
end
