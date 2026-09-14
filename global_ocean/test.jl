using TimestepperTestCases
using Oceananigans
using Oceananigans.Units
using CUDA
using Printf

# The cases of Table 1 on the eORCA025 tripolar mesh, less `WRK3-UP`, whose tracer advection this case fixes.
#
#   AB2-SE, WRK3-SE-SM05, WRK3-SE, SRK3-SE, MRK4-SE, WRK3-IM
#
# Same machinery and the same empirical time steps as `near_global/test.jl`, on the global grid.

arch = GPU()

# One grid per case: the vertical coordinate is a `MutableVerticalDiscretization`, so `ηⁿ`, the four `σ`
# scalings and `∂t_σ` live on the grid itself, and a shared one hands each case the surface state the previous
# one ended on. `global_ocean` and `global_ocean_cost` build their own grid per call.

#####
##### Phase 1: cost. Every case for 200 untimed warm-up steps and 1000 timed steps, with the substeps of its
##### nominal time step, no output and no dissipation budget, written to `global_ocean_cost_<label>_cost.jld2`.
#####

# only the timings are kept, so that each simulation leaves the GPU before the next one is built
costs = map(near_global_discretizations()) do d
    result = global_ocean_cost(d; arch)
    cost = (; result.label, Δt = result.nominal_timestep, result.seconds_per_step, result.simulated_years_per_day)
    result = nothing
    GC.gc(true)
    return cost
end

# SYPD at the production time step of each case, not at the reduced step the cost runs take
@info "[global_ocean] cost summary"
for c in costs
    @info @sprintf("  %-14s  Δt = %-12s %8.4f s/step  %8.3f SYPD", c.label, prettytime(c.Δt), c.seconds_per_step,
                   c.simulated_years_per_day)
end
flush(stderr)

#####
##### Phase 2: physics. 60 days of cold start, the cost segment to 365 days and the buoyancy-variance budget to
##### 720 days, under the `global_ocean` prefix.
#####

for d in near_global_discretizations()
    global_ocean(d; arch)
    GC.gc(true)
end
