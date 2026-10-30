using Oceananigans.BuoyancyFormulations
using Oceananigans.Models

"""
    compute_dissipation!(dissipation, model)

Compute the numerical dissipation of buoyancy variance from the cached advective fluxes, at the end of each step.

$(SIGNATURES)

The dissipation is computed using the formulation:

    A = 2 * δb★ * F - U δb²  # For advective dissipation

where `F` is the advective buoyancy flux, `U` the transport that advects the tracers, `b★ = (bⁿ⁺¹ + bⁿ) / 2` and
`b² = bⁿ⁺¹ bⁿ`. `F` and `U` are those that advanced the tracers over the step: for AB2,
`F = (1.5 + χ) Fⁿ - (0.5 + χ) Fⁿ⁻¹` and `U = (1.5 + χ) Uⁿ - (0.5 + χ) Uⁿ⁻¹`; for the Runge-Kutta schemes the flux and
transport of the final stage, or their quadrature over the stages.

This method is adapted from the exact variance budget methodology described in the paper's appendix.
"""
function compute_dissipation!(dissipation, model)
    timestepper = model.timestepper
    end_of_step(timestepper, model.clock.stage) || return nothing
    getadvection(model.advection, :T) isa Nothing && return nothing

    grid = model.grid
    C₁, C₂ = flux_combination(timestepper, eltype(grid))

    P     = dissipation.advective_production
    state = dissipation.previous_state
    Fⁿ    = dissipation.advective_fluxes.Fⁿ
    Fⁿ⁻¹  = dissipation.advective_fluxes.Fⁿ⁻¹
    bⁿ⁺¹  = Models.buoyancy_operation(model)

    launch!(architecture(grid), grid, :xyz, _assemble_advective_dissipation!,
            P, grid, C₁, C₂, Fⁿ, Fⁿ⁻¹, state.Uⁿ, state.Uⁿ⁻¹, bⁿ⁺¹, state.cⁿ⁻¹)

    return nothing
end

flux_combination(ts::QuasiAdamsBashforth2TimeStepper, FT) = (convert(FT, 1.5 + ts.χ), convert(FT, 0.5 + ts.χ))
flux_combination(::MultiStageTimeStepper, FT) = (one(FT), zero(FT))
