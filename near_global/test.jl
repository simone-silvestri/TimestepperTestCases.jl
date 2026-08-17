using TimestepperTestCases
using Oceananigans
using Oceananigans.Units
using CUDA

# The cases of Table 1 at 1/4°, defined once in `discretizations()` and shared by every test case.
#
#   WRK3-SE, WRK3-SE-SM05, WRK3-IM, SRK3-SE, AB2-SE, WRK3-UP, MRK4-SE
#
# The time steps here are empirical, not derived: on a global domain both c₁ and Δx vary, so the binding c₁ k
# is a maximum over the globe rather than a single number. The three-stage value is measured and the others
# are scaled from it by the ratio of the imaginary-axis limits, so the ratios match the idealized cases.

arch = GPU()
grid = TimestepperTestCases.near_global_grid(arch)

for d in discretizations()
    sim = TimestepperTestCases.near_global(d; arch, grid)
end

# Cost comparison across the same discretizations, each at its own time step (Section 5.4.1).
# results = run_near_global_cost(; arch, stop_time = 365days)
