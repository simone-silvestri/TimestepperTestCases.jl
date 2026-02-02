using Oceananigans.Models.VarianceDissipationComputations: assemble_advective_dissipation!
using Oceananigans.BuoyancyFormulations

"""
    compute_dissipation!(dissipation, model)

Compute the numerical dissipation of buoyancy variance from previously calculated advective fluxes.

$(SIGNATURES)

# Arguments
- `dissipation`: `BuoyancyVarianceDissipation` object containing cached fluxes
- `model`: The model containing tracers and timestepper information

# Returns
- `nothing` (modifies `dissipation.advective_production` in place)

The dissipation is computed using the formulation:

    A = 2 * δc★ * F - U δc²  # For advective dissipation

where `F` is the flux associated with advection, `U` is the advecting velocity,
`c★` and `c²` are functions of the buoyancy field.

For multi-stage timesteppers, the fluxes and velocities must account for the substepping
procedure. For AB2: `F = 1.5 Fⁿ - 0.5 Fⁿ⁻¹` and `U = 1.5 Uⁿ - 0.5 Uⁿ⁻¹`.
For RK schemes, the fluxes from the appropriate substep are used.

This method is adapted from the exact variance budget methodology described in the paper's appendix.
"""
function compute_dissipation!(dissipation, model)

    grid = model.grid

    # General velocities
    Uⁿ   = dissipation.previous_state.Uⁿ
    Uⁿ⁻¹ = dissipation.previous_state.Uⁿ⁻¹

    cⁿ⁺¹ = BuoyancyFormulations.buoyancy(model.buoyancy, model.grid, model.tracers)
    cⁿ   = dissipation.previous_state.cⁿ⁻¹

    substep = model.clock.stage
    scheme  = getadvection(model.advection, :T)

    ####
    #### Assemble the advective dissipation
    ####

    P    = dissipation.advective_production
    Fⁿ   = dissipation.advective_fluxes.Fⁿ
    Fⁿ⁻¹ = dissipation.advective_fluxes.Fⁿ⁻¹

    !(scheme isa Nothing) &&
        assemble_advective_dissipation!(P, grid, model.timestepper, substep, Fⁿ, Fⁿ⁻¹, Uⁿ, Uⁿ⁻¹, cⁿ⁺¹, cⁿ)

    return nothing
end
