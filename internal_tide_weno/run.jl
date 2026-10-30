# The internal tide cases of `discretizations()` on the Oceananigans worktree that merges `ss/sm05-barotropic` and
# `ss/weno-ghost-stencil` (see Project.toml). NumericalEarth does not support Oceananigans 0.116, so only the sources
# the internal tide needs are loaded, without `using TimestepperTestCases`.
#
#   julia --project=internal_tide_weno internal_tide_weno/run.jl            # every case, in sequence
#   julia --project=internal_tide_weno internal_tide_weno/run.jl WRK3-SE    # one case

module TimestepperTestCases

using DocStringExtensions
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.Models
using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces:
    WideTrigAveragingKernel, OptimizedAsymmetricAveragingKernel,
    LowDissipationAveragingKernel,
    ForwardBackwardScheme, RungeKutta3Scheme, averaging_shape_function,
    FrozenSlowForcing, StageQuadraticSlowForcing, ProgressiveSlowForcing
using Oceananigans.TimeSteppers: ModifiedRungeKutta4TimeStepper
using KernelAbstractions: @kernel, @index
using Printf

import Oceananigans: prognostic_state, restore_prognostic_state!
using Oceananigans.TimeSteppers: SSPRungeKuttaTimeStepper

prognostic_state(::SSPRungeKuttaTimeStepper) = nothing
restore_prognostic_state!(restored::SSPRungeKuttaTimeStepper, ::Nothing) = restored

const wall_clock = Ref(time_ns())

function print_progress(sim)
    u, v, w = sim.model.velocities
    progress = 100 * (time(sim) / sim.stop_time)
    elapsed = (time_ns() - wall_clock[]) / 1e9

    @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
            progress, iteration(sim), prettytime(sim), prettytime(elapsed),
            maximum(abs, u), maximum(abs, v), maximum(abs, w), prettytime(sim.Δt))

    wall_clock[] = time_ns()

    return nothing
end

const tracer_advection = WENO(order=7)

const source = joinpath(@__DIR__, "..", "src")

include(joinpath(source, "scheme_stability.jl"))
include(joinpath(source, "advection_schemes.jl"))
include(joinpath(source, "discretizations.jl"))
include(joinpath(source, "internal_tide.jl"))

end # module

using .TimestepperTestCases: discretizations, discretization, internal_tide, internal_tide_grid

cd(@__DIR__)

cases = isempty(ARGS) ? discretizations() : discretization.(ARGS)

for d in cases
    @info "Running the internal tide with $(d.label)"
    internal_tide(d; grid = internal_tide_grid())
end
