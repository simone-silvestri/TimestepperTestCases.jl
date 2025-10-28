using Oceananigans: fields
using Oceananigans.Grids: topology, Flat
using Oceananigans.BoundaryConditions

# Store advective and diffusive fluxes for dissipation computation
function cache_fluxes!(dissipation, model, tracer_name)
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
    cache_fluxes!(dissipation, model)

    return nothing
end

function flux_parameters(grid)
    Nx, Ny, Nz = size(grid)
    TX, TY, TZ = topology(grid)
    Fx = ifelse(TX == Flat, 1:1, 1:Nx+1)
    Fy = ifelse(TY == Flat, 1:1, 1:Ny+1)
    Fz = ifelse(TZ == Flat, 1:1, 1:Nz+1)
    return KernelParameters(Fx, Fy, Fz)
end

function cache_fluxes!(dissipation, model)

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
    advection = getadvection(model.advection, tracer_name)

    cache_advective_fluxes!(Fⁿ, Fⁿ⁻¹, grid, params, timestepper, stage, advection, U, model.buoyancy, C)


    cⁿ⁺¹ = Oceananigans.BuoyancyFormulations.buoyancy(model.buoyancy, model.grid, model.tracers)

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