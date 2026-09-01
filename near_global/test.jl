using TimestepperTestCases
using Oceananigans
using Oceananigans.Units
using CUDA

# The cases of Table 1 at 1/4°, less `WRK3-UP`, whose tracer advection this case fixes -- see
# `near_global_discretizations`.
#
#   WRK3-SE, WRK3-SE-SM05, WRK3-IM, SRK3-SE, AB2-SE, MRK4-SE
#
# The time steps here are empirical, not derived: on a global domain both c₁ and Δx vary, so the binding c₁ k
# is a maximum over the globe rather than a single number. The three-stage value is measured and the others
# are scaled from it by the ratio of the imaginary-axis limits, so the ratios match the idealized cases.

arch = GPU()

# One grid per case. The vertical coordinate is a `MutableVerticalDiscretization`, so `ηⁿ`, the four `σ`
# scalings and `∂t_σ` live on the grid itself and are written in place by whichever model holds it: sharing one
# grid across the loop hands each case the surface state the previous one ended on, and a case that ends in a
# NaN poisons every case after it at iteration zero. `regrid_bathymetry` caches the bottom height, so building
# the grid again per case costs a file read.
for d in near_global_discretizations()
    grid = TimestepperTestCases.near_global_grid(arch)
    sim = TimestepperTestCases.near_global(d; arch, grid)
end

# Cost comparison across the same discretizations, each at its own time step (Section 5.4.1).
# results = run_near_global_cost(; arch, stop_time = 365days)
