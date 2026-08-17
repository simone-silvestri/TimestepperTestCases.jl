using TimestepperTestCases
using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces

# grid
grid = TimestepperTestCases.internal_tide_grid()

# free surfaces
rk3 = SplitExplicitFreeSurfaces.RungeKutta3Scheme()

# Slow-forcing treatment across the barotropic sub-cycle. The quadratic through the three baroclinic RK
# stage values does not change the order, but removes the barotropic--baroclinic resonance near
# μ₀ = k c₀ Δt ≈ 5.7 at which the frozen forcing loses about 40% of its damping.
quad = SplitExplicitFreeSurfaces.StageQuadraticSlowForcing()

fsI = ImplicitFreeSurface()
fs1 = SplitExplicitFreeSurface(grid, substeps = 64, averaging_kernel = SplitExplicitFreeSurfaces.LowDissipationAveragingKernel())
fs2 = SplitExplicitFreeSurface(grid, substeps = 64, averaging_kernel = SplitExplicitFreeSurfaces.OptimizedAsymmetricAveragingKernel())
fs3 = SplitExplicitFreeSurface(grid, substeps = 64, timestepper = rk3, averaging_kernel = SplitExplicitFreeSurfaces.LowDissipationAveragingKernel())
fs4 = SplitExplicitFreeSurface(grid, substeps = 64, timestepper = rk3, averaging_kernel = SplitExplicitFreeSurfaces.OptimizedAsymmetricAveragingKernel())
fs5 = SplitExplicitFreeSurface(grid, substeps = 64, timestepper = rk3, averaging_kernel = SplitExplicitFreeSurfaces.OptimizedAsymmetricAveragingKernel())

# Same four, with the slow forcing reconstructed instead of frozen.
fs1q = SplitExplicitFreeSurface(grid, substeps = 64, averaging_kernel = SplitExplicitFreeSurfaces.LowDissipationAveragingKernel(), slow_forcing = quad)
fs2q = SplitExplicitFreeSurface(grid, substeps = 64, averaging_kernel = SplitExplicitFreeSurfaces.OptimizedAsymmetricAveragingKernel(), slow_forcing = quad)
fs3q = SplitExplicitFreeSurface(grid, substeps = 64, timestepper = rk3, averaging_kernel = SplitExplicitFreeSurfaces.LowDissipationAveragingKernel(), slow_forcing = quad)
fs4q = SplitExplicitFreeSurface(grid, substeps = 64, timestepper = rk3, averaging_kernel = SplitExplicitFreeSurfaces.OptimizedAsymmetricAveragingKernel(), slow_forcing = quad)

# sim = TimestepperTestCases.internal_tide(:QuasiAdamsBashforth2, free_surface = fs1)
# sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fsI, free_surface_name = "implicit")
# sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs1, free_surface_name = "split_explicit_fs1")
# sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs2, free_surface_name = "split_explicit_fs2")
# sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs3, free_surface_name = "split_explicit_fs3")
# sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs4, free_surface_name = "split_explicit_fs4")

# # Quadratic slow forcing: each of the four paired with its frozen-forcing counterpart above.
# sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs1q, free_surface_name = "split_explicit_fs1_quad")
# sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs2q, free_surface_name = "split_explicit_fs2_quad")
# sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs3q, free_surface_name = "split_explicit_fs3_quad")
# sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs4q, free_surface_name = "split_explicit_fs4_quad")
sim = TimestepperTestCases.internal_tide(:SSPRungeKutta3, free_surface = fs5, free_surface_name = "split_explicit_fs5")

# RK3-UP: featured RK3-SE config (fs4) but with 3rd-order upwind tracer advection (diffusive-spatial reference)
# sim = TimestepperTestCases.internal_tide(:SplitRungeKutta3, free_surface = fs4, free_surface_name = "up3", tracer_advection = UpwindBiased(order = 3))
