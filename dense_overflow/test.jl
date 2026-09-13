using TimestepperTestCases
using Oceananigans

# The cases of Table 1, defined once in `discretizations()` and shared by every test case.
#
#   AB2-SE, WRK3-SE-SM05, WRK3-SE, SRK3-SE, MRK4-SE, WRK3-IM, WRK3-UP

# 200 x 100 over 200 km by 2000 m, so the whole sweep runs on a CPU in minutes.
arch = CPU()

for d in discretizations()
    sim = dense_overflow(d; arch)
end
