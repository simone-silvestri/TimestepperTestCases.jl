using TimestepperTestCases
using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces

# grid
grid = TimestepperTestCases.internal_seiche_grid()

# barotropic substeppers
fb  = SplitExplicitFreeSurfaces.ForwardBackwardScheme()
rk3 = SplitExplicitFreeSurfaces.RungeKutta3Scheme()

# slow-forcing treatment across the barotropic sub-cycle
quad = SplitExplicitFreeSurfaces.StageQuadraticSlowForcing()

fs(substep, slow_forcing) =
    SplitExplicitFreeSurface(grid, substeps = 48, timestepper = substep,
                             averaging_kernel = SplitExplicitFreeSurfaces.OptimizedAsymmetricAveragingKernel(),
                             slow_forcing = slow_forcing)

# The two effects live in different ranges of μ₀ = k c₀ Δt, so both sweeps are needed.
#
#   order       : small μ₀, error against the exact modal solution. Forward--backward gives global
#                 order 1, the three-stage substep gives 2. The slow forcing is irrelevant here.
#   resonance   : μ₀ ≈ 5.7, where the frozen forcing loses about a third of its damping. The
#                 stage-value reconstruction removes it and damps ≈ 2.4× more.
μ₀_order     = (0.2, 0.4, 0.8, 1.6)
μ₀_resonance = (3.0, 4.0, 5.0, 5.5, 5.7, 6.0, 6.5, 7.0, 9.0, 12.0)

# --- order: forward--backward against the three-stage substep, frozen forcing throughout
for μ₀ in μ₀_order
    TimestepperTestCases.internal_seiche(:SplitRungeKutta3; μ₀, stop_iteration = 400,
                                         free_surface = fs(fb, FrozenSlowForcing()),
                                         free_surface_name = "fb_frozen")
    TimestepperTestCases.internal_seiche(:SplitRungeKutta3; μ₀, stop_iteration = 400,
                                         free_surface = fs(rk3, FrozenSlowForcing()),
                                         free_surface_name = "rk3_frozen")
end

# --- resonance: frozen against the stage-value quadratic, three-stage substep
for μ₀ in μ₀_resonance
    TimestepperTestCases.internal_seiche(:SplitRungeKutta3; μ₀, stop_iteration = 2000,
                                         free_surface = fs(rk3, FrozenSlowForcing()),
                                         free_surface_name = "rk3_frozen")
    TimestepperTestCases.internal_seiche(:SplitRungeKutta3; μ₀, stop_iteration = 2000,
                                         free_surface = fs(rk3, quad),
                                         free_surface_name = "rk3_quad")
end

# --- the same comparison with forward--backward substeps, where the reconstruction is expected to do
#     much less: the substep's own error dominates
for μ₀ in (5.0, 5.7, 6.0)
    TimestepperTestCases.internal_seiche(:SplitRungeKutta3; μ₀, stop_iteration = 2000,
                                         free_surface = fs(fb, FrozenSlowForcing()),
                                         free_surface_name = "fb_frozen")
    TimestepperTestCases.internal_seiche(:SplitRungeKutta3; μ₀, stop_iteration = 2000,
                                         free_surface = fs(fb, quad),
                                         free_surface_name = "fb_quad")
end

# ── MRK4: the four-stage composition of section 4 ────────────────────────────────────────────────────
#
# The seiche is the only case here that sweeps μ₀ directly, so it is where the section-4 claims can be
# measured rather than assumed: that the progressive reconstruction closes the weak baroclinic band the
# four-stage scheme carries below μ₀ ≈ 7, and that the frozen four-stage scheme is unstable over
# 4.9 < μ₀ < 8.9. The sweep is run at fixed μ₀ across schemes, so the comparison is at equal barotropic
# phase per step rather than at equal time step.
mrk4 = ModifiedRungeKutta4TimeStepper()
prog = ProgressiveSlowForcing(mrk4)

# The band the progressive reconstruction has to close.
μ₀_band = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)

for μ₀ in μ₀_band
    TimestepperTestCases.internal_seiche(:ModifiedRungeKutta4; μ₀, stop_iteration = 2000,
                                         free_surface = fs(rk3, prog),
                                         free_surface_name = "mrk4_progressive")
    TimestepperTestCases.internal_seiche(:SplitRungeKutta3; μ₀, stop_iteration = 2000,
                                         free_surface = fs(rk3, quad),
                                         free_surface_name = "rk3_quad")
end
