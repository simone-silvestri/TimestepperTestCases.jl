module TimestepperTestCases

export internal_tide, idealized_coast, channel_simulation

using DocStringExtensions
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.Models
using KernelAbstractions: @kernel, @index
using Printf

wall_clock = Ref(time_ns())

"""
    print_progress(sim)

Print simulation progress information including completion percentage, iteration number,
simulation time, wall clock time, maximum velocity components, and next time step.

$(SIGNATURES)

# Arguments
- `sim`: The `Simulation` object to print progress for.

# Returns
- `nothing`

This callback function is typically added to simulations using `add_callback!` to monitor
progress during long-running simulations.
"""
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

# For all test cases
const tracer_buffer_scheme = WENO(order=5, buffer_scheme=Centered())
const tracer_advection     = WENO(order=7, buffer_scheme=tracer_buffer_scheme)

include("BuoyancyVarianceDissipationComputations/BuoyancyVarianceDissipationComputations.jl")
using .BuoyancyVarianceDissipationComputations: BuoyancyVarianceDissipation

# Simulations!
include("internal_tide.jl")
include("idealized_coast.jl")
include("channel_simulation.jl")

using Oceananigans
using Oceananigans.AbstractOperations: grid_metric_operation
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.ImmersedBoundaries: column_depthᶜᶜᵃ,
                                       static_column_depthᶜᶜᵃ,
                                       column_depthᶠᶜᵃ,
                                       static_column_depthᶠᶜᵃ,
                                       column_depthᶜᶠᵃ,
                                       static_column_depthᶜᶠᵃ

using Oceananigans.Operators: volume
using Oceananigans.Utils
using Statistics

using KernelAbstractions: @kernel, @index

# Diagnostics!
include("diagnostics.jl")
include("load_idealized_coast_case.jl")
include("load_internal_tide_case.jl")
include("load_channel_case.jl")

end
