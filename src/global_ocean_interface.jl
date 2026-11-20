function global_ocean end
function global_ocean_grid end
function global_ocean_vertical_discretization end
function global_ocean_barotropic_rate end
function global_ocean_cost end

"""
$(TYPEDSIGNATURES)

The discretizations integrated for five years in the global configuration: every manuscript discretization except
`WRK3-IM`, which is only timed.
"""
global_ocean_discretizations() = filter(d -> d.label != "WRK3-IM", discretizations())

"""
$(TYPEDSIGNATURES)

Simulated years per wall-clock day, given the production timestep [s] and measured wall-clock seconds per step.
"""
global_ocean_simulated_years_per_day(nominal_timestep, seconds_per_step) = nominal_timestep / seconds_per_step / 365
