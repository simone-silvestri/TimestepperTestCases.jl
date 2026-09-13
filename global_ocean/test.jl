using TimestepperTestCases
using Oceananigans
using Oceananigans.Units
using CUDA

# The cases of Table 1 on the eORCA025 tripolar mesh, less `WRK3-UP`, whose tracer advection this case fixes.
#
#   AB2-SE, WRK3-SE-SM05, WRK3-SE, SRK3-SE, MRK4-SE, WRK3-IM
#
# Same machinery and the same empirical time steps as `near_global/test.jl`, on the global grid and under the
# `global_ocean` output prefix.

arch = GPU()

# One grid per case: the vertical coordinate is a `MutableVerticalDiscretization`, so `ηⁿ`, the four `σ`
# scalings and `∂t_σ` live on the grid itself, and a shared one hands each case the surface state the previous
# one ended on. `global_ocean` builds its own grid per call.
for d in near_global_discretizations()
    sim = global_ocean(d; arch)
end

# Cost comparison across the same discretizations (Section 5.4.1). `run_near_global_cost` runs with
# `dissipation = false`, so the timings are the model's and not the buoyancy-variance budget's.
# results = run_near_global_cost(; arch, case = global_ocean, stop_time = 365days)
