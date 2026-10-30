using TimestepperTestCases
using Oceananigans

# WRK3-SE-SM05 with the First-Stage coupling of Ducousso et al. in place of the per-stage one (see
# `src/first_stage_coupling.jl`). Writes `internal_tide_SplitRungeKutta3_WRK3-FS-SM05.jld2`.

grid = TimestepperTestCases.internal_tide_grid()
sim  = TimestepperTestCases.internal_tide(first_stage_discretization(); grid)
