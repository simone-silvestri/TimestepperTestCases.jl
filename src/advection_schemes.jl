#####
##### The WENO reconstructions shared by every case
#####
##### Next to walls, domain buffers and immersed boundaries the stencils are continued past the last active node,
##### so every reconstruction keeps its order up to the boundary.
#####

# A single `WENO`: an adaptive implicit `time_discretization` acts on the vertical flux alone
weno_tracer_advection(order, time_discretization) = WENO(; order, time_discretization)

function split_momentum_advection(order, time_discretization)

    vorticity_order, remaining_order = isnothing(order) ? (9, 5) : (order, order)

    vorticity_scheme               = WENO(order=vorticity_order)
    divergence_scheme              = WENO(order=remaining_order)
    kinetic_energy_gradient_scheme = WENO(order=remaining_order)
    vertical_advection_scheme      = WENO(; order = remaining_order, time_discretization)

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
