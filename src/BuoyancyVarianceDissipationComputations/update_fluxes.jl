using Oceananigans: fields
using Oceananigans.Grids: topology, Flat
using Oceananigans.Utils
using Oceananigans.BoundaryConditions
using Oceananigans.BuoyancyFormulations
using Oceananigans.Models

"""
    cache_fluxes!(dissipation, model)

Cache advective fluxes and update velocity fields for dissipation computation.

$(SIGNATURES)

# Arguments
- `dissipation`: `BuoyancyVarianceDissipation` object
- `model`: The model containing velocities and tracers

# Returns
- `nothing` (modifies `dissipation` in place)

This function updates the velocity fields (`Uⁿ`, `Uⁿ⁻¹`) and caches the advective buoyancy
fluxes (`Fⁿ`, `Fⁿ⁻¹`) needed for computing dissipation at the next time step. The caching
procedure depends on the timestepper type and stage.
"""
function cache_fluxes!(dissipation, model)
    grid = model.grid
    sz   = size(model.tracers[1].data)
    of   = model.tracers[1].data.offsets

    params = KernelParameters(sz, of)

    Uⁿ   = dissipation.previous_state.Uⁿ
    Uⁿ⁻¹ = dissipation.previous_state.Uⁿ⁻¹
    U    = model.velocities
    timestepper = model.timestepper
    stage = model.clock.stage

    update_transport!(Uⁿ, Uⁿ⁻¹, grid, params, timestepper, stage, U)
    finally_cache_fluxes!(dissipation, model)

    return nothing
end

"""
    flux_parameters(grid)

Compute kernel parameters for flux computation based on grid topology.

$(SIGNATURES)

# Arguments
- `grid`: Grid object

# Returns
- `KernelParameters` object with appropriate index ranges for flux computation

This function determines the correct index ranges for computing fluxes at face locations,
accounting for flat (1D/2D) dimensions in the grid topology.
"""
function flux_parameters(grid)
    Nx, Ny, Nz = size(grid)
    TX, TY, TZ = topology(grid)
    Fx = ifelse(TX == Flat, 1:1, 1:Nx+1)
    Fy = ifelse(TY == Flat, 1:1, 1:Ny+1)
    Fz = ifelse(TZ == Flat, 1:1, 1:Nz+1)
    return KernelParameters(Fx, Fy, Fz)
end

"""
    finally_cache_fluxes!(dissipation, model)

Cache advective buoyancy fluxes from temperature and salinity tracers.

$(SIGNATURES)

# Arguments
- `dissipation`: `BuoyancyVarianceDissipation` object
- `model`: The model containing tracers and advection scheme

# Returns
- `nothing` (modifies `dissipation.advective_fluxes` in place)

This function computes the advective fluxes of buoyancy by combining temperature and
salinity fluxes using the linear equation of state. The fluxes are computed at face
locations and cached for use in dissipation computation.
"""
function finally_cache_fluxes!(dissipation, model)

    # Grab tracer properties
    C    = model.tracers
    cⁿ⁻¹ = dissipation.previous_state.cⁿ⁻¹

    grid = model.grid
    arch = architecture(grid)
    U = model.velocities
    params = flux_parameters(grid)
    stage  = model.clock.stage
    timestepper = model.timestepper

    ####
    #### Update the advective fluxes and compute gradient squared
    ####

    Fⁿ   = dissipation.advective_fluxes.Fⁿ
    Fⁿ⁻¹ = dissipation.advective_fluxes.Fⁿ⁻¹
    advection = getadvection(model.advection, :T)

    cache_advective_fluxes!(Fⁿ, Fⁿ⁻¹, grid, params, timestepper, stage, advection, U, model.buoyancy, C)


    cⁿ⁺¹ = Models.buoyancy_operation(model.buoyancy, model.grid, model.tracers)

    if timestepper isa QuasiAdamsBashforth2TimeStepper
        set!(cⁿ⁻¹, cⁿ⁺¹)
        fill_halo_regions!(cⁿ⁻¹)
    elseif (timestepper isa RungeKuttaScheme) && (stage == length(timestepper.β))
        set!(cⁿ⁻¹, cⁿ⁺¹)
        fill_halo_regions!(cⁿ⁻¹)
    end

    return nothing
end

cache_advective_fluxes!(Fⁿ, Fⁿ⁻¹, grid, params, ::QuasiAdamsBashforth2TimeStepper, stage, advection, U, b, C) =
    launch!(architecture(grid), grid, params, _cache_advective_fluxes!, Fⁿ, Fⁿ⁻¹, grid, advection, U, b, C)

function cache_advective_fluxes!(Fⁿ, Fⁿ⁻¹, grid, params, ts::SplitRungeKuttaTimeStepper, stage, advection, U, b, C)
    if stage == length(ts.β)-1
        launch!(architecture(grid), grid, params, _cache_advective_fluxes!, Fⁿ, grid, advection, U, b, C)
    end
end

update_transport!(Uⁿ, Uⁿ⁻¹, grid, params, ::QuasiAdamsBashforth2TimeStepper, stage, U) =
    launch!(architecture(grid), grid, params, _update_transport!, Uⁿ, Uⁿ⁻¹, grid, U)

function update_transport!(Uⁿ, Uⁿ⁻¹, grid, params, ts::SplitRungeKuttaTimeStepper, stage, U)
    if stage == length(ts.β)-1
        launch!(architecture(grid), grid, params, _update_transport!, Uⁿ, grid, U)
    end
end

@kernel function _update_transport!(Uⁿ, Uⁿ⁻¹, grid, U)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        Uⁿ⁻¹.u[i, j, k] = Uⁿ.u[i, j, k]
        Uⁿ⁻¹.v[i, j, k] = Uⁿ.v[i, j, k]
        Uⁿ⁻¹.w[i, j, k] = Uⁿ.w[i, j, k]
          Uⁿ.u[i, j, k] = U.u[i, j, k] * Axᶠᶜᶜ(i, j, k, grid)
          Uⁿ.v[i, j, k] = U.v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid)
          Uⁿ.w[i, j, k] = U.w[i, j, k] * Azᶜᶜᶠ(i, j, k, grid)
    end
end

@kernel function _update_transport!(Uⁿ, grid, U)
    i, j, k = @index(Global, NTuple)

    @inbounds Uⁿ.u[i, j, k] = U.u[i, j, k] * Axᶠᶜᶜ(i, j, k, grid)
    @inbounds Uⁿ.v[i, j, k] = U.v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid)
    @inbounds Uⁿ.w[i, j, k] = U.w[i, j, k] * Azᶜᶜᶠ(i, j, k, grid)
end