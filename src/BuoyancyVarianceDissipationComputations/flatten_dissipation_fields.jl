"""
    flatten_dissipation_fields(t::BuoyancyVarianceDissipation)

Flatten the dissipation fields into a named tuple for output.

$(SIGNATURES)

# Arguments
- `t`: `BuoyancyVarianceDissipation` object

# Returns
- Named tuple containing `Abx`, `Aby`, `Abz` fields representing advective buoyancy dissipation
  in the x, y, and z directions respectively

This function extracts the dissipation fields from the internal C-grid vector structure
into a flat named tuple suitable for use as simulation outputs.
"""
function flatten_dissipation_fields(t::BuoyancyVarianceDissipation) 
    Abx = t.advective_production.x
    Aby = t.advective_production.y
    Abz = t.advective_production.z

    return (; Abx, Aby, Abz)
end
