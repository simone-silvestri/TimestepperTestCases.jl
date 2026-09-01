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

timesteppers = Dict("WRK3-SE" => "SplitRungeKutta3",
                    "SRK3-SE" => "SSPRungeKutta3",
                    "MRK4-SE" => "ModifiedRungeKutta4")

for factor in (1//4, 1//8), label in ("WRK3-SE", "SRK3-SE", "MRK4-SE")
    tag  = label * "-dt" * string(denominator(factor))
    file = "internal_tide/internal_tide_" * timesteppers[label] * "_" * tag * ".jld2"

    if isfile(file)
        @info "skipping $tag, output already present"
        continue
    end

    sim = TimestepperTestCases.internal_tide(discretization(label);
                                             timestep_factor = factor, label = tag)
end