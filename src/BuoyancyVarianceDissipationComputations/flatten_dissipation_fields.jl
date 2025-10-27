"""
    flatten_dissipation_fields(t::VarianceDissipation)

Flatten the dissipation fields of a `VarianceDissipation` object into a named tuple containing:

- The dissipation associated with the advection scheme in fields named `A-tracername-dir`
- The dissipation associated with the closures in fields names `D-tracername-dir`
"""
function flatten_dissipation_fields(t::BuoyancyVarianceDissipation) 
    Abx = t.advective_production.x
    Aby = t.advective_production.y
    Abz = t.advective_production.z

    return (; Abx, Aby, Abz)
end
