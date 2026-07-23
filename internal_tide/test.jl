using TimestepperTestCases
using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces

# grid
grid = TimestepperTestCases.internal_tide_grid()

# free surfaces
rk3 = SplitExplicitFreeSurfaces.RungeKutta3Scheme()

fsI = ImplicitFreeSurface()
fs1 = SplitExplicitFreeSurface(grid, substeps = 60, averaging_kernel = SplitExplicitFreeSurfaces.LowDissipationAveragingKernel())
fs2 = SplitExplicitFreeSurface(grid, substeps = 60, averaging_kernel = SplitExplicitFreeSurfaces.WideTrig74AveragingKernel())
fs3 = SplitExplicitFreeSurface(grid, substeps = 60, timestepper = rk3, averaging_kernel = SplitExplicitFreeSurfaces.LowDissipationAveragingKernel())
fs4 = SplitExplicitFreeSurface(grid, substeps = 60, timestepper = rk3, averaging_kernel = SplitExplicitFreeSurfaces.WideTrig74AveragingKernel())

sim = TimestepperTestCases.internal_tide(:QuasiAdamsBashforth2, free_surface = fs1)
sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fsI, free_surface_name = "implicit")
sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs1, free_surface_name = "split_explicit_fs1")
sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs2, free_surface_name = "split_explicit_fs2")
sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs3, free_surface_name = "split_explicit_fs3")
sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs4, free_surface_name = "split_explicit_fs4")

# RK3-UP: featured RK3-SE config (fs4) but with 3rd-order upwind tracer advection (diffusive-spatial reference)
sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs4, free_surface_name = "up3", tracer_advection = UpwindBiased(order = 3))