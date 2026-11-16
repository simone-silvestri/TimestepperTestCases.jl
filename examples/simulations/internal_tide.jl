using TimestepperTestCases
using Oceananigans
using Oceananigans.Units

# Output is saved every hour, counted in steps for AB2-SE: an output aligned in time shortens the step that lands on
# it, and QAB2 restarts a shortened step with forward Euler.
function internal_tide_output_schedule(d, factor = 1)
    d.timestepper == :QuasiAdamsBashforth2 || return TimeInterval(1hour)
    Δt = factor * TimestepperTestCases.internal_tide_timestep(Val(d.timestepper))
    return IterationInterval(round(Int, 1hour / Δt))
end

for d in discretizations()
    grid = TimestepperTestCases.internal_tide_grid()
    TimestepperTestCases.internal_tide(d; grid, output_schedule = internal_tide_output_schedule(d))
end

# ══ Temporal-resolution test (section 5.1) ════════════════════════════════════════════════════════════════════

for factor in (1//4, 1//8), label in ("AB2-SE", "WRK3-SE", "SRK3-SE", "MRK4-SE")
    case = discretization(label)
    tag  = label * "-dt" * string(denominator(factor))
    file = joinpath("internal_tide", "internal_tide_$(case.timestepper)_$(tag).jld2")

    if isfile(file)
        @info "skipping $tag, output already present"
        continue
    end

    TimestepperTestCases.internal_tide(case; timestep_factor = factor, label = tag,
                                       output_schedule = internal_tide_output_schedule(case, factor))
end
