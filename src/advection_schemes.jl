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

    vorticity_scheme               = WENO(order=vorticity_order, minimum_buffer_upwind_order=1)
    divergence_scheme              = WENO(order=remaining_order, minimum_buffer_upwind_order=1)
    kinetic_energy_gradient_scheme = WENO(order=remaining_order, minimum_buffer_upwind_order=1)
    vertical_advection_scheme      = WENO(; order = remaining_order, time_discretization, stretched=(:z,),
                                            minimum_buffer_upwind_order=1)

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
