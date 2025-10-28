using Oceananigans.Models.VarianceDissipationComputations: assemble_advective_dissipation!
using Oceananigans.BuoyancyFormulations

"""
    compute_dissipation!(dissipation, model)

Compute the numerical dissipation for tracer `tracer_name`, from the previously calculated advective and diffusive fluxes,
the formulation is:

    A = 2 * δc★ * F - U δc²  # For advective dissipation
    D = 2 * δc★ * F          # For diffusive dissipation

Where ``F'' is the flux associated with the particular process,``U'' is the advecting velocity,
while ``c★'' and ``c²'' are functions defined above.
Note that ``F'' and ``U'' need to be numerically accurate for the budgets to close,
i.e. for the AB2 scheme:

    F = 1.5 Fⁿ - 0.5 Fⁿ⁻¹
    U = 1.5 Uⁿ - 0.5 Uⁿ⁻¹

For a RK3 method (not implemented at the moment), the whole substepping procedure needs to be accounted for.
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
