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
#####   :ghost_cells  `GhostCells()`, which keeps the full order and completes the stencil with ghost values that
#####                 blend the mirror image of the active run with its quadratic extrapolation
#####   :upwind       first-order upwind, monotone, in exactly those cells and nowhere else
#####   :default      `Centered(order=2)` for tracers and `UpwindBiased(order=1)` for momentum, the Oceananigans
#####                 defaults
#####

using Oceananigans.Advection: GhostCells
using Oceananigans.Advection.Adapt: Adapt
using Oceananigans.Architectures: on_architecture

# Oceananigans does not adapt `GhostCells`, so the `RefValue` Δt of an adaptive implicit `scheme` reaches GPU kernels
Adapt.adapt_structure(to, scheme::GhostCells) =
    GhostCells(Adapt.adapt(to, scheme.scheme), scheme.curvature_weight, scheme.monotone)

Oceananigans.Architectures.on_architecture(to, scheme::GhostCells) =
    GhostCells(on_architecture(to, scheme.scheme), scheme.curvature_weight, scheme.monotone)

"""
    const default_tracer_boundary_scheme = :ghost_cells

The reconstruction every case terminates its tracer WENO chains in unless the caller says otherwise.
"""
const default_tracer_boundary_scheme = :ghost_cells

"""
    const default_momentum_boundary_scheme = :upwind

The reconstruction every case terminates its momentum WENO chains in unless the caller says otherwise.
"""
const default_momentum_boundary_scheme = :upwind

"""
    tracer_boundary_reconstruction(boundary_scheme)
    momentum_boundary_reconstruction(boundary_scheme)

The reconstruction in which a tracer, respectively a momentum, WENO chain terminates. `boundary_scheme` is a
`Val` of `:default`, `:upwind` or `:ghost_cells`, as returned by [`boundary_scheme_value`](@ref); `nothing`
leaves the chain on the Oceananigans default.
"""
tracer_boundary_reconstruction(::Val{:default})     = nothing
tracer_boundary_reconstruction(::Val{:upwind})      = UpwindBiased(order=1)
tracer_boundary_reconstruction(::Val{:ghost_cells}) = GhostCells()

momentum_boundary_reconstruction(::Val{:default})     = UpwindBiased(order=1)
momentum_boundary_reconstruction(::Val{:upwind})      = UpwindBiased(order=1)
momentum_boundary_reconstruction(::Val{:ghost_cells}) = GhostCells()

"""
    boundary_scheme_value(boundary_scheme)

`Val(boundary_scheme)`, after checking that it names one of the three reconstructions. A `Val` passes through,
so a case function may be handed either the symbol or the value it dispatches on.
"""
function boundary_scheme_value(boundary_scheme::Symbol)
    boundary_scheme ∈ (:default, :upwind, :ghost_cells) ||
        throw(ArgumentError("boundary_scheme must be :default, :upwind or :ghost_cells, got $boundary_scheme"))

    return Val(boundary_scheme)
end

boundary_scheme_value(boundary_scheme::Val) = boundary_scheme

"""
    boundary_scheme_name(boundary_scheme)

The symbol a `Val` boundary scheme carries: `:default`, `:upwind` or `:ghost_cells`.
"""
boundary_scheme_name(::Val{name}) where name = name

"""
    boundary_scheme_suffix(tracer_boundary_scheme, momentum_boundary_scheme;
                           reference_tracer_scheme = default_tracer_boundary_scheme,
                           reference_momentum_scheme = default_momentum_boundary_scheme)

Filename tag naming the boundary reconstructions, empty at `reference_tracer_scheme` and
`reference_momentum_scheme`, the defaults of the case, so that a run at the default keeps the plain filename. The two are
named separately only where they disagree, one scheme applied to both being the common case.

Two of the three options reach the model as the same advection type and differ only in the reconstruction the
chain terminates in, so a filename built from the type of the scheme cannot tell them apart.
"""
function boundary_scheme_suffix(tracer_boundary_scheme, momentum_boundary_scheme;
                                reference_tracer_scheme = default_tracer_boundary_scheme,
                                reference_momentum_scheme = default_momentum_boundary_scheme)

    tracer   = boundary_scheme_name(boundary_scheme_value(tracer_boundary_scheme))
    momentum = boundary_scheme_name(boundary_scheme_value(momentum_boundary_scheme))

    tracer === reference_tracer_scheme && momentum === reference_momentum_scheme && return ""
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
