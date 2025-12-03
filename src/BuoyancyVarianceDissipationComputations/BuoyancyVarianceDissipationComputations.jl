module BuoyancyVarianceDissipationComputations

export BuoyancyVarianceDissipation

using Oceananigans.Advection
using Oceananigans.BoundaryConditions
using Oceananigans.Grids: architecture, AbstractGrid
using Oceananigans.Utils
using Oceananigans.Fields
using Oceananigans.Fields: Field, VelocityFields
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.TimeSteppers: QuasiAdamsBashforth2TimeStepper,
                                 SplitRungeKuttaTimeStepper

using Oceananigans.TurbulenceClosures: viscosity,
                                       diffusivity,
                                       ScalarDiffusivity,
                                       ScalarBiharmonicDiffusivity,
                                       AbstractTurbulenceClosure,
                                       HorizontalFormulation,
                                       _diffusive_flux_x,
                                       _diffusive_flux_y,
                                       _diffusive_flux_z

using Oceananigans.Advection: _advective_tracer_flux_x,
                              _advective_tracer_flux_y,
                              _advective_tracer_flux_z

using Oceananigans: UpdateStateCallsite
using Oceananigans.Operators: volume
using Oceananigans.Utils: IterationInterval, ConsecutiveIterations
using KernelAbstractions: @kernel, @index

const RungeKuttaScheme = SplitRungeKuttaTimeStepper

struct BuoyancyVarianceDissipation{P, A, S}
    advective_production :: P
    advective_fluxes :: A
    previous_state :: S
end

"""
    c_grid_vector(grid)

Create a named tuple of C-grid vector fields (x, y, z components at face locations).

$(SIGNATURES)

# Arguments
- `grid`: Grid on which to create the fields

# Returns
- Named tuple with `x`, `y`, `z` fields at XFace, YFace, and ZFace locations respectively

This helper function creates the vector fields needed to store advective fluxes and
dissipation rates at their proper staggered grid locations.
"""
function c_grid_vector(grid)
    x = XFaceField(grid)
    y = YFaceField(grid)
    z = ZFaceField(grid)
    return (; x, y, z)
end

"""
    BuoyancyVarianceDissipation(grid; Uⁿ⁻¹, Uⁿ)

Construct a `BuoyancyVarianceDissipation` object for computing buoyancy variance dissipation diagnostics.

$(SIGNATURES)

# Arguments
- `grid`: The grid on which the buoyancy field is defined

# Keyword Arguments
- `Uⁿ⁻¹`: The velocity field at the previous time step (default: `VelocityFields(grid)`)
- `Uⁿ`: The velocity field at the current time step (default: `VelocityFields(grid)`)

# Returns
- `BuoyancyVarianceDissipation` object ready to be used as a simulation callback

This function computes the variance dissipation diagnostics for buoyancy, which is computed
from temperature and salinity tracers using the linear equation of state. The diagnostics
include the numerical dissipation implicit to the advection scheme.

This diagnostic is especially useful for models that use a dissipative advection scheme
like WENO or UpwindBiased, as it quantifies the spurious mixing introduced by numerical
discretization.

!!! compat "Time stepper compatibility"
    The variance dissipation diagnostic is supported for `QuasiAdamsBashforth2TimeStepper`
    and `SplitRungeKuttaTimeStepper`.
"""
function BuoyancyVarianceDissipation(grid; 
                                     Uⁿ⁻¹ = VelocityFields(grid),
                                     Uⁿ   = VelocityFields(grid))

    P    = c_grid_vector(grid)
    Fⁿ   = c_grid_vector(grid)
    Fⁿ⁻¹ = c_grid_vector(grid)
    cⁿ⁻¹ = CenterField(grid)
    
    previous_state   = (; cⁿ⁻¹, Uⁿ⁻¹, Uⁿ)
    advective_fluxes = (; Fⁿ, Fⁿ⁻¹)

    return BuoyancyVarianceDissipation(P, advective_fluxes, previous_state)
end

"""
    (ϵ::BuoyancyVarianceDissipation)(model)

Compute buoyancy variance dissipation and update flux caches for the next time step.

$(SIGNATURES)

# Arguments
- `model`: The `HydrostaticFreeSurfaceModel` containing temperature and salinity tracers

# Returns
- `nothing` (modifies `ϵ` in place)

This function is called as a simulation callback to compute the advective dissipation
of buoyancy variance at each time step. It first computes dissipation from previously
cached fluxes, then updates the flux caches for the next time step.

The model must have `:T` and `:S` tracers to compute buoyancy.
"""
function (ϵ::BuoyancyVarianceDissipation)(model)

    # Check if the model has a velocity field
    if !hasproperty(model, :velocities)
        throw(ArgumentError("Model must have a velocity field."))
    end

    # Check if the model has tracers
    if !(hasproperty(model.tracers, :T) && hasproperty(model.tracers, :S))
        throw(ArgumentError("Model must have tracers T and S."))
    end

    # First we compute the dissipation from previously computed fluxes
    compute_dissipation!(ϵ, model)

    # Then we update the fluxes to be used in the next time step
    cache_fluxes!(ϵ, model)

    return nothing
end

@inline getadvection(advection, tracer_name) = advection
@inline getadvection(advection::NamedTuple, tracer_name) = @inbounds advection[tracer_name]

const f = Face()
const c = Center()

include("update_fluxes.jl")
include("advective_dissipation.jl")
include("compute_dissipation.jl")
include("flatten_dissipation_fields.jl")

import Oceananigans.Simulations: Callback, validate_schedule

# Specific `Callback` for `VarianceDissipation` computations.
# A VarianceDissipation object requires a `ConsecutiveIteration` schedule to make sure
# that the computed fluxes are correctly used in the next time step.
# Also, the `VarianceDissipation` object needs to be called on `UpdateStateStepCallsite` to be correct.
function Callback(func::BuoyancyVarianceDissipation, schedule=IterationInterval(1);
                  parameters = nothing,
                  callsite = UpdateStateCallsite())

    if !(callsite isa UpdateStateCallsite)
        @warn "BuoyancyVarianceDissipation callback must be called on UpdateStateCallsite. Changing `callsite` to `UpdateStateCallsite()`."
        callsite = UpdateStateCallsite()
    end

    schedule = validate_schedule(func, schedule)

    return Callback(func, schedule, callsite, parameters)
end

validate_schedule(::BuoyancyVarianceDissipation, schedule) = throw(ArgumentError("the provided schedule $schedule is not supported for VarianceDissipation computations. \n" *
                                                                         "Use an `IterationInterval` schedule instead."))

function validate_schedule(::BuoyancyVarianceDissipation, schedule::IterationInterval)
    if !(schedule == IterationInterval(1))
        @warn "BuoyancyVarianceDissipation callback must be called every Iteration or on `ConsecutiveIterations`. \n" *
              "Changing `schedule` to `ConsecutiveIterations(schedule)`."
        schedule = ConsecutiveIterations(schedule)
    end
    return schedule
end

function validate_schedule(::BuoyancyVarianceDissipation, schedule::ConsecutiveIterations)
    if !(schedule.parent isa IterationInterval)
       throw(ArgumentError("the provided schedule $schedule is not supported for BuoyancyVarianceDissipation computations. \n" *
                                                                         "Use an `IterationInterval` schedule instead."))
    end
    return schedule
end

end # module
