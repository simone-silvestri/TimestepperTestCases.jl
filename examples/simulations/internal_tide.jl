using TimestepperTestCases
using Oceananigans

for d in discretizations()
    grid = TimestepperTestCases.internal_tide_grid()
    TimestepperTestCases.internal_tide(d; grid)
end

for factor in (1//4, 1//8), label in ("WRK3-SE", "SRK3-SE", "MRK4-SE")
    case = discretization(label)
    tag  = label * "-dt" * string(denominator(factor))
    file = joinpath("internal_tide", "internal_tide_$(case.timestepper)_$(tag).jld2")

    if isfile(file)
        @info "skipping $tag, output already present"
        continue
    end

    TimestepperTestCases.internal_tide(case; timestep_factor = factor, label = tag)
end
