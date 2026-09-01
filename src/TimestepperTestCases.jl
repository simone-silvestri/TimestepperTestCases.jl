module TimestepperTestCases

export internal_tide, internal_seiche, idealized_coast, channel_simulation
export near_global, near_global_grid, near_global_discretizations, run_near_global_cost
export load_near_global, load_near_global_cases, near_global_cost_table
export near_global_diffusivity_profile, near_global_diffusivity_map, near_global_diffusivity_hovmoller
export near_global_surface_speed, near_global_surface_kinetic_energy, near_global_eddy_kinetic_energy
export near_global_zonal_spectrum

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

export stability_limit, baroclinic_timestep, barotropic_substeps, barotropic_courant
export Discretization, discretizations, discretization, timestep_and_free_surface, timestep_ratio, plot_style
export ModifiedRungeKutta4TimeStepper, ProgressiveSlowForcing

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

# Upstream `prognostic_state(::AbstractTimeStepper)` reads `Gⁿ` and `G⁻`, fields the SSP stepper does not
# carry; like `SplitRungeKuttaTimeStepper`, it holds no state between time steps, so checkpointing skips it.
import Oceananigans: prognostic_state, restore_prognostic_state!
using Oceananigans.TimeSteppers: SSPRungeKuttaTimeStepper

prognostic_state(::SSPRungeKuttaTimeStepper) = nothing
restore_prognostic_state!(restored::SSPRungeKuttaTimeStepper, ::Nothing) = restored

# For all test cases
const tracer_buffer_scheme = WENO(order=5, buffer_scheme=Centered())
const tracer_advection     = WENO(order=7, buffer_scheme=tracer_buffer_scheme)

include("scheme_stability.jl")
include("discretizations.jl")

include("BuoyancyVarianceDissipationComputations/BuoyancyVarianceDissipationComputations.jl")
using .BuoyancyVarianceDissipationComputations: BuoyancyVarianceDissipation

# Simulations!
include("internal_tide.jl")
include("internal_seiche.jl")
include("idealized_coast.jl")
include("channel_simulation.jl")

# Realistic near-global quarter-degree cost + buoyancy-dissipation case (needs NumericalEarth).
include("near_global.jl")

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
include("load_internal_seiche_case.jl")
include("load_channel_case.jl")
include("load_near_global_case.jl")

end
