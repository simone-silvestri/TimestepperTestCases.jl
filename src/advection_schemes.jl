#####
##### The WENO reconstructions shared by every case
#####
##### Next to walls, domain buffers and immersed boundaries the tracer stencils are continued past the last
##### active node. The momentum reconstructions instead step down to first-order upwind at the boundary.
#####

# A single `WENO`: an adaptive implicit `time_discretization` acts on the vertical flux alone
weno_tracer_advection(order, time_discretization) = WENO(; order, time_discretization, stretched=(:z,))

function split_momentum_advection(order, time_discretization)

    vorticity_order, remaining_order = isnothing(order) ? (9, 5) : (order, order)

    # Momentum only: against a wall the stencil steps down through its own family and terminates in
    # first-order upwind, which is what `WENOVectorInvariant` applies and what the earlier campaign ran. The
    # tracer reconstruction keeps its own buffer and is deliberately not changed here.
    momentum_buffer = Oceananigans.Advection.DecreasingOrderAdvectionScheme(
                          boundary_scheme = Oceananigans.Advection.UpwindBiased(order = 1))

    vorticity_scheme               = WENO(order=vorticity_order, buffer_scheme=momentum_buffer)
    divergence_scheme              = WENO(order=remaining_order, buffer_scheme=momentum_buffer)
    kinetic_energy_gradient_scheme = WENO(order=remaining_order, buffer_scheme=momentum_buffer)
    vertical_advection_scheme      = WENO(; order = remaining_order, time_discretization, stretched=(:z,),
                                            buffer_scheme=momentum_buffer)

    return VectorInvariant(; vorticity_scheme,
                             vertical_advection_scheme,
                             divergence_scheme,
                             kinetic_energy_gradient_scheme)
end

function split_flux_form_momentum_advection(order, time_discretization)
    horizontal = WENO(; order)
    vertical   = WENO(; order, time_discretization)
    return FluxFormAdvection(horizontal, horizontal, vertical)
end
