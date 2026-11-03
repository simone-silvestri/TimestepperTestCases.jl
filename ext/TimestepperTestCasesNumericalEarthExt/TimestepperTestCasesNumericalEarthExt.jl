module TimestepperTestCasesNumericalEarthExt

using DocStringExtensions
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.Models
using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces:
    WideTrigAveragingKernel, OptimizedAsymmetricAveragingKernel,
    LowDissipationAveragingKernel,
    ForwardBackwardScheme, RungeKutta3Scheme, averaging_shape_function,
    FrozenSlowForcing, ProgressiveSlowForcing
using Oceananigans.TimeSteppers: ModifiedRungeKutta4TimeStepper
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Utils
using KernelAbstractions: @kernel, @index
using Statistics
using Printf

using TimestepperTestCases
using TimestepperTestCases: BuoyancyVarianceDissipation, Discretization, baroclinic_timestep, barotropic_speed, barotropic_substeps,
                            campaign_timestepper, forwarded_tracer_advection, split_momentum_advection, staggered_wavenumber, weno_tracer_advection,
                            global_ocean_discretizations, global_ocean_simulated_years_per_day

import TimestepperTestCases: global_ocean, global_ocean_grid, global_ocean_vertical_discretization, global_ocean_barotropic_rate,
                             global_ocean_cost

include("global_simulation.jl")
include("global_ocean.jl")

end # module
