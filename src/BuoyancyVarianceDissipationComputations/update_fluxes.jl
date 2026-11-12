using Oceananigans: fields
using Oceananigans.Grids: topology, Flat
using Oceananigans.Utils
using Oceananigans.BoundaryConditions
using Oceananigans.BuoyancyFormulations
using Oceananigans.Models

"""
    cache_fluxes!(dissipation, model)

Cache the advective buoyancy flux and the transport of the stage that has just ended.

$(SIGNATURES)

A `HydrostaticFreeSurfaceModel` builds the transport that advects its tracers inside each stage, after the
barotropic step, from the velocities the previous stage ended with. The flux of a stage is therefore computed at its
end, from that transport and from the tracers and the grid stretching the previous stage ended with. AB2 combines the
flux of the step with the one of the previous step, the split Runge-Kutta schemes keep the flux of the final stage,
and the strong-stability-preserving scheme the quadrature `Σₘ βₘ Fᵐ` over the stages.
"""
function cache_fluxes!(dissipation, model)
    grid        = model.grid
    timestepper = model.timestepper
    weights     = stage_weights(timestepper, model.clock.stage, eltype(grid))
    isnothing(weights) && return nothing

    weight, keep = weights
    store_previous_advective_fluxes!(dissipation, timestepper)

    state  = dissipation.previous_state
    params = flux_parameters(grid)
    U      = model.transport_velocities

    launch!(architecture(grid), grid, params, _accumulate_transport!, state.Uⁿ, grid, U, state.σᵐ⁻¹, weight, keep)

    Fⁿ           = dissipation.advective_fluxes.Fⁿ
    tracers      = (T = state.Tᵐ⁻¹, S = state.Sᵐ⁻¹)
    coefficients = buoyancy_coefficients(model.buoyancy)

    for name in (:T, :S)
        advection = getadvection(model.advection, name)
        launch!(architecture(grid), grid, params, _accumulate_advective_fluxes!,
                Fⁿ, grid, advection, U, tracers[name], state.σᵐ⁻¹, weight * coefficients[name], keep)
        keep = one(keep)
    end

    return nothing
end

"""
    stage_weights(timestepper, stage, FT)

The weight with which the flux of `stage` enters the flux of the step, and whether the sum already under way is
kept, or `nothing` if that flux does not enter the step.
"""
stage_weights(::QuasiAdamsBashforth2TimeStepper, stage, FT) = (one(FT), zero(FT))

stage_weights(ts::SplitRungeKuttaTimeStepper, stage, FT) = stage == ts.Nstages ? (one(FT), zero(FT)) : nothing

function stage_weights(ts::SSPRungeKuttaTimeStepper, stage, FT)
    β = ssp_quadrature_weights(ts.coefficients)
    return convert(FT, β[stage]), convert(FT, stage == 1 ? 0 : 1)
end

end_of_step(::QuasiAdamsBashforth2TimeStepper, stage) = true
end_of_step(ts::MultiStageTimeStepper, stage) = stage == ts.Nstages

store_previous_advective_fluxes!(dissipation, timestepper) = nothing

function store_previous_advective_fluxes!(dissipation, ::QuasiAdamsBashforth2TimeStepper)
    Fⁿ, Fⁿ⁻¹ = dissipation.advective_fluxes
    Uⁿ, Uⁿ⁻¹ = dissipation.previous_state.Uⁿ, dissipation.previous_state.Uⁿ⁻¹

    for (previous, current) in zip((Fⁿ⁻¹..., Uⁿ⁻¹...), (Fⁿ..., Uⁿ...))
        parent(previous) .= parent(current)
    end

    return nothing
end

"""
    store_stage_state!(dissipation, model)

Keep the tracers and the horizontal grid stretching the next stage starts from and, at the end of the step, the
buoyancy the next step starts from.
"""
function store_stage_state!(dissipation, model)
    grid  = model.grid
    state = dissipation.previous_state

    parent(state.Tᵐ⁻¹) .= parent(model.tracers.T)
    parent(state.Sᵐ⁻¹) .= parent(model.tracers.S)
    launch!(architecture(grid), grid, horizontal_flux_parameters(grid), _store_horizontal_stretching!, state.σᵐ⁻¹, grid)

    if end_of_step(model.timestepper, model.clock.stage)
        set!(state.cⁿ⁻¹, Models.buoyancy_operation(model))
        fill_halo_regions!(state.cⁿ⁻¹)
    end

    return nothing
end

function face_ranges(grid)
    Nx, Ny, Nz = size(grid)
    TX, TY, TZ = topology(grid)
    Fx = ifelse(TX == Flat, 1:1, 1:Nx+1)
    Fy = ifelse(TY == Flat, 1:1, 1:Ny+1)
    Fz = ifelse(TZ == Flat, 1:1, 1:Nz+1)
    return Fx, Fy, Fz
end

"""
    flux_parameters(grid)

Kernel parameters spanning every face of `grid`, with flat directions reduced to a single index.
"""
flux_parameters(grid) = KernelParameters(face_ranges(grid)...)

function horizontal_flux_parameters(grid)
    Fx, Fy, _ = face_ranges(grid)
    return KernelParameters(Fx, Fy)
end
