#####
##### The reconstruction in which a WENO buffer chain terminates
#####
##### Every case reconstructs with a WENO chain, and there is always one cell in which the stencil no longer
##### fits: the domain buffer and, on an `ImmersedBoundaryGrid`, any cell adjacent to an `inactive_node`
##### (`Advection/immersed_advective_fluxes.jl`). The interior keeps the full order, so what happens in that one
##### cell is a choice of its own. It is made once here and shared by every case.
#####
##### The three options are
#####
#####   :cwenoz   the third-order central-WENO reconstruction of Semplice, Travaglia and Puppo (2022), whose
#####             stencil extends only inwards, blending an inward parabola, a linear polynomial and a constant
#####             with Z-weights
#####   :upwind   first-order upwind, monotone, in exactly those cells and nowhere else
#####   :default  `Centered(order=2)` for tracers and `UpwindBiased(order=1)` for momentum, the Oceananigans
#####             defaults
#####
##### The tracers terminate in `:cwenoz` and the momentum in `:upwind`, the two being chosen separately because
##### the CWENOZ scales are well posed for one and not for the other. A tracer is a single field with a single
##### gradient per direction, so `reference_gradient` carries its units and the stencil estimate of `ϵ` tracks
##### the column. The horizontal vector-invariant terms instead reconstruct a vorticity, a divergence flux and
##### a squared velocity, three different units, so no reference gradient is available to them and `d⁰` sits
##### pinned at its cap; on the near-global bathymetry that termination loses the zonal velocity within the
##### first day, whereas the same reconstruction on the tracers alone runs.
#####

using Oceananigans.Advection: CWENOZ

"""
    const default_tracer_boundary_scheme = :cwenoz

The reconstruction every case terminates its tracer WENO chains in unless the caller says otherwise.
"""
const default_tracer_boundary_scheme = :cwenoz

"""
    const default_momentum_boundary_scheme = :upwind

The reconstruction every case terminates its momentum WENO chains in unless the caller says otherwise.
"""
const default_momentum_boundary_scheme = :upwind

"""
    tracer_boundary_reconstruction(boundary_scheme, reference_gradient)
    momentum_boundary_reconstruction(boundary_scheme, reference_gradient)

The reconstruction in which a tracer, respectively a momentum, WENO chain terminates. `boundary_scheme` is a
`Val` of `:default`, `:upwind` or `:cwenoz`, as returned by [`boundary_scheme_value`](@ref); `nothing` leaves
the chain on the Oceananigans default.

`reference_gradient` sets the oscillation scale `ϵ = (∇ref Δ)²` below which CWENOZ reads the data as smooth and
recovers third order. It carries the units of the reconstructed field per unit length, so it belongs to one
variable and one direction; zero estimates it from the stencil as `min(I¹, I¹')`, the smaller of the two linear
oscillations, which is the local increment squared.

The stencil estimate sets `ϵ` to the local oscillation itself, so at a genuine step `ϵ` grows with the roughness,
`τ/ϵ` stays near one and the constant candidate never takes over: the blend keeps third order across the
discontinuity and undershoots. The constant is capped at `d⁰ = 0.01` against `d° = 0.74`, so it wins only once
`τ = 5/3 c²` reaches ~74 `ϵ`, that is once the second difference `c` exceeds ~6.7 `∇ref Δ`. With `∇ref` at the
typical gradient of the field, `∇ref Δ` is the typical *first* difference and a genuine step gets 2% constant
weight -- no limiting at all, plain third order. Firing the constant needs `∇ref` about a decade below the
typical gradient, which is where the horizontal tracer scales of the near-global case sit, small enough that a
topographic step activates the constant while smooth data is untouched.

The vertical directions keep the stencil estimate, because `ϵ ∝ Δ²` grows monotonically downwards while the
per-cell increment `∇ψ Δz` does not, so a single `∇ref` sets scales differing by orders of magnitude down one
column. The horizontal momentum terms reconstruct three different units, so they carry no reference gradient
either, and there the estimate is the only option available.
"""
tracer_boundary_reconstruction(::Val{:default}, reference_gradient) = nothing
tracer_boundary_reconstruction(::Val{:upwind},  reference_gradient) = UpwindBiased(order=1)
tracer_boundary_reconstruction(::Val{:cwenoz},  reference_gradient) = CWENOZ(; reference_gradient)

momentum_boundary_reconstruction(::Val{:default}, reference_gradient) = UpwindBiased(order=1)
momentum_boundary_reconstruction(::Val{:upwind},  reference_gradient) = UpwindBiased(order=1)
momentum_boundary_reconstruction(::Val{:cwenoz},  reference_gradient) = CWENOZ(; reference_gradient)

"""
    boundary_scheme_value(boundary_scheme)

`Val(boundary_scheme)`, after checking that it names one of the three reconstructions. A `Val` passes through,
so a case function may be handed either the symbol or the value it dispatches on.
"""
function boundary_scheme_value(boundary_scheme::Symbol)
    boundary_scheme ∈ (:default, :upwind, :cwenoz) ||
        throw(ArgumentError("boundary_scheme must be :default, :upwind or :cwenoz, got $boundary_scheme"))

    return Val(boundary_scheme)
end

boundary_scheme_value(boundary_scheme::Val) = boundary_scheme

"""
    boundary_scheme_name(boundary_scheme)

The symbol a `Val` boundary scheme carries: `:default`, `:upwind` or `:cwenoz`.
"""
boundary_scheme_name(::Val{name}) where name = name

"""
    boundary_scheme_suffix(tracer_boundary_scheme, momentum_boundary_scheme)

Filename tag naming the boundary reconstructions, empty at [`default_tracer_boundary_scheme`](@ref) and
[`default_momentum_boundary_scheme`](@ref) so that a run at the default keeps the plain filename. The two are
named separately only where they disagree, one scheme applied to both being the common case.

Two of the three options reach the model as the same advection type and differ only in the reconstruction the
chain terminates in, so a filename built from the type of the scheme cannot tell them apart.
"""
function boundary_scheme_suffix(tracer_boundary_scheme, momentum_boundary_scheme)
    tracer   = boundary_scheme_name(boundary_scheme_value(tracer_boundary_scheme))
    momentum = boundary_scheme_name(boundary_scheme_value(momentum_boundary_scheme))

    tracer === default_tracer_boundary_scheme && momentum === default_momentum_boundary_scheme && return ""
    tracer === momentum && return "_$tracer"

    return "_tracer_$(tracer)_momentum_$(momentum)"
end

"""
    split_flux_form_advection(boundary_reconstruction, order, time_discretization, boundary_scheme,
                              horizontal_reference_gradient, vertical_reference_gradient)

Flux-form advection of order `order` whose horizontal and vertical reconstructions terminate in the boundary
schemes `boundary_reconstruction` returns for `horizontal_reference_gradient` and `vertical_reference_gradient`
respectively. Only the vertical direction takes `time_discretization`, which is where the vertically implicit
treatment applies.
"""
function split_flux_form_advection(boundary_reconstruction, order, time_discretization, boundary_scheme,
                                   horizontal_reference_gradient, vertical_reference_gradient)

    horizontal_boundary_scheme = boundary_reconstruction(boundary_scheme, horizontal_reference_gradient)
    vertical_boundary_scheme   = boundary_reconstruction(boundary_scheme, vertical_reference_gradient)

    horizontal = WENO(; order, boundary_scheme = horizontal_boundary_scheme)
    vertical   = WENO(; order, time_discretization, boundary_scheme = vertical_boundary_scheme)

    return FluxFormAdvection(horizontal, horizontal, vertical)
end

"""
    split_tracer_advection(order, time_discretization, boundary_scheme,
                           horizontal_reference_gradient, vertical_reference_gradient)

Tracer advection of order `order` whose reconstructions terminate in the tracer boundary schemes carrying the
two reference gradients, one per direction.
"""
split_tracer_advection(order, time_discretization, boundary_scheme,
                       horizontal_reference_gradient, vertical_reference_gradient) =
    split_flux_form_advection(tracer_boundary_reconstruction, order, time_discretization, boundary_scheme,
                              horizontal_reference_gradient, vertical_reference_gradient)

"""
    split_momentum_advection(order, time_discretization, boundary_scheme,
                             horizontal_reference_gradient, vertical_reference_gradient)

Vector-invariant momentum advection whose four reconstructions terminate in a boundary scheme chosen per
direction, reproducing `WENOVectorInvariant` in every other respect -- `order = nothing` keeps its per-term
defaults, a ninth-order vorticity flux and fifth order everywhere else.

The vorticity, divergence and kinetic-energy-gradient terms are the horizontal ones, and they reconstruct a
vorticity, a divergence flux and a squared velocity: three different units, so `horizontal_reference_gradient`
is dimensionally meaningful only at zero, where the scale is read off the stencil. The vertical term
reconstructs velocity, so `vertical_reference_gradient` is a shear in inverse seconds. Both default to zero for
the reason given in [`tracer_boundary_reconstruction`](@ref).
"""
function split_momentum_advection(order, time_discretization, boundary_scheme,
                                  horizontal_reference_gradient, vertical_reference_gradient)

    vorticity_order, remaining_order = isnothing(order) ? (9, 5) : (order, order)

    horizontal_boundary_scheme = momentum_boundary_reconstruction(boundary_scheme, horizontal_reference_gradient)
    vertical_boundary_scheme   = momentum_boundary_reconstruction(boundary_scheme, vertical_reference_gradient)

    vorticity_scheme               = WENO(order=vorticity_order, boundary_scheme=horizontal_boundary_scheme)
    divergence_scheme              = WENO(order=remaining_order, boundary_scheme=horizontal_boundary_scheme)
    kinetic_energy_gradient_scheme = WENO(order=remaining_order, boundary_scheme=horizontal_boundary_scheme)
    vertical_advection_scheme      = WENO(; order = remaining_order, time_discretization,
                                            boundary_scheme = vertical_boundary_scheme)

    return VectorInvariant(; vorticity_scheme,
                             vertical_advection_scheme,
                             divergence_scheme,
                             kinetic_energy_gradient_scheme)
end

"""
    split_flux_form_momentum_advection(order, time_discretization, boundary_scheme,
                                       horizontal_reference_gradient, vertical_reference_gradient)

Flux-form momentum advection, for the internal tide, which advects momentum in flux form rather than in the
vector-invariant form the channel and the near-global case use. The reconstructions are those of
[`split_tracer_advection`](@ref) built on the momentum boundary reconstruction, so `:default` terminates in a
first-order upwind rather than in a second-order centered.

Both reference gradients are shears in inverse seconds, momentum being what is reconstructed in every
direction, contrarily to the vector-invariant form.
"""
split_flux_form_momentum_advection(order, time_discretization, boundary_scheme,
                                   horizontal_reference_gradient, vertical_reference_gradient) =
    split_flux_form_advection(momentum_boundary_reconstruction, order, time_discretization, boundary_scheme,
                              horizontal_reference_gradient, vertical_reference_gradient)

"""
    tracer_advection_scheme(boundary_scheme, default; order, time_discretization,
                            horizontal_reference_gradient, vertical_reference_gradient)

`default` where `boundary_scheme` is `:default`, and [`split_tracer_advection`](@ref) otherwise.

The case's own scheme is threaded through rather than rebuilt, so that `:default` reproduces it bit for bit:
the shared `TimestepperTestCases.tracer_advection` terminates its WENO7 chain at order five, one step earlier
than the chain `WENO(; order, boundary_scheme)` builds, and the idealized cases have output on disk that was
produced with it.
"""
tracer_advection_scheme(::Val{:default}, default; kw...) = default

tracer_advection_scheme(boundary_scheme::Val, default;
                        order = 7,
                        time_discretization = ExplicitTimeDiscretization(),
                        horizontal_reference_gradient = 0,
                        vertical_reference_gradient = 0) =
    split_tracer_advection(order, time_discretization, boundary_scheme,
                           horizontal_reference_gradient, vertical_reference_gradient)
