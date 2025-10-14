module TimestepperTestCases

export internal_tide, baroclinic_adjustment, overflow

using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.Models
using Printf

wall_clock = Ref(time_ns())

function print_progress(sim)
    model = sim.model
    u, v, w = model.velocities
    progress = 100 * (time(sim) / sim.stop_time)
    elapsed = (time_ns() - TimestepperTestCases.wall_clock[]) / 1e9

    @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
            progress, iteration(sim), prettytime(sim), prettytime(elapsed),
            maximum(abs, u), maximum(abs, v), maximum(abs, w), prettytime(sim.Δt))

    TimestepperTestCases.wall_clock[] = time_ns()

    return nothing
end

include("baroclinic_adjustment.jl")
include("internal_tide.jl")
include("idealized_coast.jl")

end
