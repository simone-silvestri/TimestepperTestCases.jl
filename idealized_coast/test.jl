using TimestepperTestCases
using Oceananigans
using CUDA

# The cases of Table 1, plus the SM05 reference and the four-stage composition of section 4, defined once in
# `discretizations()` and shared by every test case.
#
#   WRK3-SE, WRK3-SE-SM05, WRK3-IM, SRK3-SE, AB2-SE, WRK3-UP, MRK4-SE

arch = GPU()

# `lowres = true` is the configuration the runs use: 96 x 96 over 192 km, i.e. the 2 km horizontal resolution
# the manuscript quotes. It also sets the time step, through k = pi/dx.
lowres = true

for d in discretizations()
    sim = idealized_coast(d; arch, lowres)
end
