#####
##### The discretizations of Table 1, defined once and shared by every test case
#####
##### A test case is the product of two independent choices: the physical configuration (internal tide,
##### channel, global) and the numerical discretization. Only the discretizations are
##### collected here; the configuration enters through its `stability_parameters`, from which the time step and
##### the substep count are *derived* rather than stored. Nothing in this file carries a time step, a substep
##### count or a ratio between them: those follow from the composition's own coefficients through
##### [`baroclinic_timestep`](@ref) and [`barotropic_substeps`](@ref), so they cannot drift from the scheme
##### they describe.
#####

"""
    struct Discretization

One row of Table 1: a baroclinic composition, a free-surface formulation, and — when the free surface is
split-explicit — the barotropic substep integrator, averaging kernel and slow-forcing treatment that go with
it, plus the tracer advection scheme.

Fields
======
- `label`: the case name, used for output filenames and for the labels in the figures
- `timestepper`: the baroclinic composition, as the symbol `HydrostaticFreeSurfaceModel` takes
- `barotropic_timestepper`: the substep integrator, or `nothing` for an implicit free surface
- `averaging_kernel`: the barotropic averaging filter, or `nothing` for an implicit free surface
- `slow_forcing`: how the slow forcing is represented across the sub-cycle
- `tracer_advection`: the tracer advection scheme
- `implicit_free_surface`: whether the free surface is solved implicitly

The time step is not a field: it follows from the composition through [`baroclinic_timestep`](@ref), which
places every scheme at the same fraction of its own stability limit -- see [`timestep_ratio`](@ref).
"""
struct Discretization{TS, BT, AK, SF, TA}
    label :: String
    timestepper :: TS
    barotropic_timestepper :: BT
    averaging_kernel :: AK
    slow_forcing :: SF
    tracer_advection :: TA
    implicit_free_surface :: Bool
end

function Discretization(label, timestepper;
                        barotropic_timestepper = RungeKutta3Scheme(),
                        averaging_kernel       = OptimizedAsymmetricAveragingKernel(),
                        slow_forcing           = FrozenSlowForcing(),
                        tracer_advection       = TimestepperTestCases.tracer_advection,
                        implicit_free_surface  = false)

    return Discretization(label, timestepper, barotropic_timestepper, averaging_kernel,
                          slow_forcing, tracer_advection, implicit_free_surface)
end

"""
    timestep_ratio(d::Discretization)

The time step of `d` relative to the WRK3 reference: the ratio of the imaginary-axis limits,
`θ★(scheme)/θ★(WRK3)`, which is 0.29 for AB2 and 1.609 for MRK4. Derived rather than stored, so it cannot
disagree with the composition it describes.
"""
timestep_ratio(d::Discretization) = stability_limit(d.timestepper) / stability_limit(:SplitRungeKutta3)

"""
    forwarded_tracer_advection(d::Discretization)

The tracer advection of `d` where it departs from the shared `TimestepperTestCases.tracer_advection`, and
`nothing` where it does not.

A `Discretization` that names no scheme of its own must not override the one of the case: only `WRK3-UP`, which
exists to carry a third-order upwind, reaches the model with the scheme it stores.
"""
forwarded_tracer_advection(d::Discretization) =
    d.tracer_advection === TimestepperTestCases.tracer_advection ? nothing : d.tracer_advection

Base.summary(d::Discretization) = string("Discretization(\"", d.label, "\")")
Base.show(io::IO, d::Discretization) = print(io, summary(d))

"""
    discretizations()

Return the seven campaign configurations. AB2-SE and WRK3-SE-SM05 use forward-backward substeps with the
Shchepetkin–McWilliams filter. WRK3-SE, SRK3-SE, MRK4-SE, and WRK3-UP use RK3 substeps with the optimized
averaging filter; WRK3-IM uses an implicit free surface. WRK3-SE, MRK4-SE, and WRK3-UP reconstruct the slow
forcing progressively. WRK3-UP uses third-order upwind tracer advection.
"""
function discretizations()
    fb   = ForwardBackwardScheme()
    sm05 = averaging_shape_function

    # Ordered as the ladder of Table 3, which changes one ingredient at a time: the AB2 incumbent, then the
    # outer-stepper swap at fixed barotropic treatment, then the barotropic treatment at fixed outer stepper.
    return [Discretization("AB2-SE", :QuasiAdamsBashforth2;
                           barotropic_timestepper = fb, averaging_kernel = sm05),

            Discretization("WRK3-SE-SM05", :SplitRungeKutta3;
                           barotropic_timestepper = fb, averaging_kernel = sm05),

            Discretization("WRK3-SE", :SplitRungeKutta3;
                           slow_forcing = ProgressiveSlowForcing()),

            Discretization("SRK3-SE", :SSPRungeKutta3),

            Discretization("MRK4-SE", :ModifiedRungeKutta4;
                           slow_forcing = ProgressiveSlowForcing(campaign_timestepper(:ModifiedRungeKutta4))),

            Discretization("WRK3-IM", :SplitRungeKutta3;
                           barotropic_timestepper = nothing, averaging_kernel = nothing,
                           implicit_free_surface  = true),

            Discretization("WRK3-UP", :SplitRungeKutta3;
                           slow_forcing     = ProgressiveSlowForcing(),
                           tracer_advection = UpwindBiased(order = 3))]
end

"""
    discretization(label)

The [`Discretization`](@ref) of `discretizations()` with the given label.
"""
function discretization(label)
    cases = discretizations()
    index = findfirst(d -> d.label == label, cases)
    isnothing(index) && throw(ArgumentError("Unknown discretization $label."))
    return cases[index]
end

"""
    timestep_and_free_surface(d::Discretization, grid, stability_parameters)

The baroclinic time step and the materialized free surface for discretization `d` on a configuration described
by `stability_parameters`, a named tuple `(; N², H, Δx)`.

Both are derived: the time step from the theoretical limit of `d.timestepper` -- on the first baroclinic mode,
or on advection where `stability_parameters` supplies a speed `U` and advection is the tighter of the two --
and the substep count as the smallest that keeps the barotropic Courant number at the `barotropic_cfl` of the substep integrator, 0.8 for forward-backward and 1.3 for the
three-stage Runge-Kutta substep ([`barotropic_cfl`](@ref)). The substep count therefore depends on the substep integrator as well as on the
time step, which is why the forward-backward cases carry more substeps than the Runge-Kutta ones.
"""
function timestep_and_free_surface(d::Discretization, grid, stability_parameters)
    Δt = baroclinic_timestep(d.timestepper; stability_parameters...)

    d.implicit_free_surface && return Δt, ImplicitFreeSurface()

    substeps = barotropic_substeps(d.barotropic_timestepper;
                                   stability_parameters.H, stability_parameters.Δx, Δt,
                                   d.averaging_kernel)

    free_surface = SplitExplicitFreeSurface(grid; substeps,
                                            averaging_kernel = d.averaging_kernel,
                                            timestepper      = d.barotropic_timestepper,
                                            slow_forcing     = d.slow_forcing)

    return Δt, free_surface
end

#####
##### Plotting style
#####

"""
    plot_style(d)   /   plot_style(label)

Colour and line style for a case, so that every Results notebook draws the same case the same way.

The palette is Okabe--Ito, which survives greyscale reproduction and the common forms of colour blindness.
Style carries meaning rather than being decorative: dashed marks the FB/SM05 barotropic treatment, dotted the
implicit free surface. WRK3-SE and SRK3-SE deliberately share a colour and differ only in dash, because they
share a stability polynomial and every other setting -- whether the two curves lie on top of each other is
itself the result, and giving them separate colours hides it.
"""
function plot_style(label::AbstractString; distinct_colors = false)
    # Where the panel already spends line style on something else -- the channel figures use it to separate two
    # operators within a case -- style is not available to tell WRK3-SE from SRK3-SE, so they need their own
    # colours. `distinct_colors = true` gives every case a colour of its own, at the cost of no longer showing
    # the agreement of the two three-stage schemes as an overlay.
    srk3_color = distinct_colors ? "#56B4E9" : "#0072B2"

    label == "AB2-SE"       && return (color = :black,      linestyle = :solid)
    label == "WRK3-SE-SM05" && return (color = "#E69F00",   linestyle = :dash)
    label == "WRK3-SE"      && return (color = "#0072B2",   linestyle = :solid)
    label == "SRK3-SE"      && return (color = srk3_color,  linestyle = :dashdot)
    label == "MRK4-SE"      && return (color = "#CC79A7",   linestyle = :dash)
    label == "WRK3-IM"      && return (color = "#D55E00",   linestyle = :dot)
    label == "WRK3-UP"      && return (color = "#009E73",   linestyle = :solid)

    throw(ArgumentError("No plot style for $label. Known: " *
                        string([d.label for d in discretizations()])))
end

plot_style(d::Discretization; kw...) = plot_style(d.label; kw...)
