using TimestepperTestCases
using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces

SplitExplicitFreeSurfaces.variable_density_barotropic_mode[] = true

# barotropic sub-step schemes
fb  = SplitExplicitFreeSurfaces.ForwardBackwardScheme()
rk3 = SplitExplicitFreeSurfaces.RungeKutta3Scheme()

# averaging filters
mu2    = SplitExplicitFreeSurfaces.LowDissipationAveragingKernel()   # μ₂ = 0
trig74 = SplitExplicitFreeSurfaces.WideTrig74AveragingKernel()       # μ₂ = μ₃ = 0 (wide window)

# same cases as internal_tide/test.jl (channel names its outputs through `testcase`):
sim = channel_simulation(timestepper = :QuasiAdamsBashforth2, testcase = "ab2")
sim = channel_simulation(timestepper = :SplitRungeKutta3, free_surface = ImplicitFreeSurface(), testcase = "implicit")
sim = channel_simulation(timestepper = :SplitRungeKutta3, barotropic_timestepper = fb,  averaging_kernel = mu2,    testcase = "split_explicit_fs1")
sim = channel_simulation(timestepper = :SplitRungeKutta3, barotropic_timestepper = fb,  averaging_kernel = trig74, testcase = "split_explicit_fs2")
sim = channel_simulation(timestepper = :SplitRungeKutta3, barotropic_timestepper = rk3, averaging_kernel = mu2,    testcase = "split_explicit_fs3")
sim = channel_simulation(timestepper = :SplitRungeKutta3, barotropic_timestepper = rk3, averaging_kernel = trig74, testcase = "split_explicit_fs4")
