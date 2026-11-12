using TimestepperTestCases
using NumericalEarth
using WorldOceanAtlasTools
using Oceananigans
using Oceananigans.Units
using CUDA
using Printf

arch = GPU()

# Retain timing records and release each simulation before constructing the next case.
costs = map(discretizations()) do d
    result = global_ocean_cost(d; arch)
    cost = (; result.label, Δt = result.nominal_timestep, result.seconds_per_step, result.simulated_years_per_day)
    result = nothing
    GC.gc(true)
    return cost
end

@info "[global_ocean] cost summary"
for c in costs
    @info @sprintf("  %-14s  Δt = %-12s %8.4f s/step  %8.3f SYPD", c.label, prettytime(c.Δt), c.seconds_per_step,
                   c.simulated_years_per_day)
end
flush(stderr)

for d in global_ocean_discretizations()
    global_ocean(d; arch)
    GC.gc(true)
end
