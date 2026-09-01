using TimestepperTestCases
using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces

# The cases of Table 1, plus the SM05 reference and the four-stage composition of section 4, defined once in
# `discretizations()` and shared by every test case.
#
#   WRK3-SE, WRK3-SE-SM05, WRK3-IM, SRK3-SE, AB2-SE, WRK3-UP, MRK4-SE
#
# The channel's stratification is surface-intensified, so its c₁ comes from the WKB integral -- see
# `channel_stability_parameters`.

for d in discretizations()[4:end]
    sim = channel_simulation(d)
end
