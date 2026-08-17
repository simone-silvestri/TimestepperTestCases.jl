using TimestepperTestCases
using Oceananigans

# The cases of Table 1, plus the SM05 reference and the four-stage composition of section 4. They are defined
# once in `discretizations()` and shared by every test case; see there for what each one isolates.
#
#   WRK3-SE, WRK3-SE-SM05, WRK3-IM, SRK3-SE, AB2-SE, WRK3-UP, MRK4-SE
#
# The time step and the barotropic substep count are derived from the discretization and from
# `internal_tide_stability_parameters()`, not prescribed here.

grid = TimestepperTestCases.internal_tide_grid()

for d in discretizations()
    sim = TimestepperTestCases.internal_tide(d; grid)
end
