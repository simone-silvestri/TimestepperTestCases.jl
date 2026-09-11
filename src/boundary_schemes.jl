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
    tracer_boundary_reconstruction(boundary_scheme)
    momentum_boundary_reconstruction(boundary_scheme)

The reconstruction in which a tracer, respectively a momentum, WENO chain terminates. `boundary_scheme` is a
`Val` of `:default`, `:upwind` or `:cwenoz`, as returned by [`boundary_scheme_value`](@ref); `nothing` leaves
the chain on the Oceananigans default.

CWENOZ reads the oscillation scale `ϵ` off the stencil, as `min(I¹, I¹')`, the smaller of the two linear
oscillations, so `ϵ` is the local increment squared. At a genuine step `ϵ` grows with the roughness, `τ/ϵ`
stays near one and the constant candidate never takes over: the blend keeps third order across the
discontinuity and undershoots.
"""
tracer_boundary_reconstruction(::Val{:default}) = nothing
tracer_boundary_reconstruction(::Val{:upwind})  = UpwindBiased(order=1)
tracer_boundary_reconstruction(::Val{:cwenoz})  = CWENOZ()

momentum_boundary_reconstruction(::Val{:default}) = UpwindBiased(order=1)
momentum_boundary_reconstruction(::Val{:upwind})  = UpwindBiased(order=1)
momentum_boundary_reconstruction(::Val{:cwenoz})  = CWENOZ()

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

function split_flux_form_advection(boundary_reconstruction, order, time_discretization, boundary_scheme)

    horizontal_boundary_scheme = boundary_reconstruction(boundary_scheme)
    vertical_boundary_scheme   = boundary_reconstruction(boundary_scheme)

    horizontal = WENO(; order, boundary_scheme = horizontal_boundary_scheme)
    vertical   = WENO(; order, time_discretization, boundary_scheme = vertical_boundary_scheme)

    return FluxFormAdvection(horizontal, horizontal, vertical)
end

split_tracer_advection(order, time_discretization, boundary_scheme) =
    split_flux_form_advection(tracer_boundary_reconstruction, order, time_discretization, boundary_scheme)

function split_momentum_advection(order, time_discretization, boundary_scheme)

    vorticity_order, remaining_order = isnothing(order) ? (9, 5) : (order, order)

    horizontal_boundary_scheme = momentum_boundary_reconstruction(boundary_scheme)
    vertical_boundary_scheme   = momentum_boundary_reconstruction(boundary_scheme)

    vorticity_scheme               = WENO(order=vorticity_order, boundary_scheme=horizontal_boundary_scheme)
    divergence_scheme              = WENO(order=remaining_order, boundary_scheme=horizontal_boundary_scheme)
    kinetic_energy_gradient_scheme = WENO(order=remaining_order, boundary_scheme=horizontal_boundary_scheme)
    vertical_advection_scheme      = WENO(; order = remaining_order, time_discretization, boundary_scheme = vertical_boundary_scheme)

    return VectorInvariant(; vorticity_scheme,
                             vertical_advection_scheme,
                             divergence_scheme,
                             kinetic_energy_gradient_scheme)
end

split_flux_form_momentum_advection(order, time_discretization, boundary_scheme) =
    split_flux_form_advection(momentum_boundary_reconstruction, order, time_discretization, boundary_scheme)

tracer_advection_scheme(::Val{:default}, default; kw...) = default

tracer_advection_scheme(boundary_scheme::Val, default;
                        order = 7,
                        time_discretization = ExplicitTimeDiscretization()) =
    split_tracer_advection(order, time_discretization, boundary_scheme)
