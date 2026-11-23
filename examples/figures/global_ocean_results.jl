# # Global ocean: numerical diffusivity
#
# In this example, we diagnose the numerical diffusivity of the six five-year global runs from their buoyancy-variance
# budgets, the `global_ocean_<label>_dissipation.jld2` files of the record. The latitude–depth section is Figure 16 of
# the manuscript. Run it with `julia --project=examples examples/figures/global_ocean_results.jl /path/to/global_ocean`.

using TimestepperTestCases
using Oceananigans
using Statistics
include("../shared/data_sources.jl") #src
include("../shared/compact_data.jl") #src
include("../shared/plotting.jl") #src
include("../shared/global_plots.jl") #src

directory = isempty(ARGS) ? joinpath(@__DIR__, "..", "..", "global_ocean") : first(ARGS)
labels = ["AB2-SE", "WRK3-SE-SM05", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-UP"]

# ## Reducing the budgets
#
# Each file holds the three-dimensional variance tendencies and squared buoyancy gradients, averaged over 90-day windows.
# We first reduce every window to sums by depth, by distance from topography and by latitude band, keeping the tendency
# and the gradient separate so that the diffusivity can be formed after summing. The staggered terms are moved to cell
# centers, and cells with nonfinite metrics are left out of both sums. The reductions are cached in `data/global_ocean/`;
# passing that directory redraws the figures without the three-dimensional fields.

budgets = map(labels) do label
    @info "Reading global buoyancy-variance budget" label
    filename = "global_ocean_$(label)_dissipation.jld2"
    compact_dataset("global_ocean", filename, prepare_global_budget; directory, reuse=true)
end

# ## Numerical diffusivity
#
# Over any region, the numerical diffusivity is
#
# ```math
# \kappa = - \frac{\sum P}{2 \sum G},
# ```
#
# where ``P`` is the variance tendency and ``G`` the squared buoyancy gradient, both already multiplied by the cell
# volume. The budget joins at day 1095, so its first window is incomplete; we use the eight complete windows, which end at
# days 1170 to 1800, in every case.

@assert all(budget -> budget.window_ends ≈ collect(1170:90:1800), budgets)
numerical_diffusivity(P, G) = @. -P / (2G)
window_sum(values) = dropdims(sum(values; dims=3); dims=3)

profiles = map(labels, budgets) do label, budget
    P = window_sum(budget.wall_tendency)
    G = window_sum(budget.wall_gradient)
    return (; label, budget.z, κ=numerical_diffusivity(vec(sum(P; dims=1)), vec(sum(G; dims=1))),
              wall_diffusivity=numerical_diffusivity(P, G))
end

plot_global_diffusivity_profile(profiles) #src
# ![Depth profiles](global_diffusivity_profile.png)
#
# Unlike in the idealized cases, the Runge–Kutta runs stay close to AB2-SE at every depth. Only WRK3-UP departs, with a
# global diffusivity about twice as large: in this eddy-permitting configuration the spatial discretization, rather
# than the time stepping, appears to control the numerical mixing.
#
# Near land the WENO stencils lose order, so we also split the profiles by distance from topography. `wall_distance`
# counts the cells to the nearest land along the grid lines, the only directions a WENO stencil reaches, from one to
# five cells and beyond.

plot_global_diffusivity_wall_distance(profiles) #src
# ![Profiles by distance from topography](global_diffusivity_wall_distance.png)

# ## Latitude–depth sections
#
# Summing over quarter-degree latitude bands instead gives the sections of Figure 16. The first panel shows WRK3-SE on a
# logarithmic scale from 10⁻⁶ to 10⁻³ m² s⁻¹, and the others the difference of each configuration from WRK3-SE.

sections = map(labels, budgets) do label, budget
    latitudes = (budget.latitude_edges[1:end-1] + budget.latitude_edges[2:end]) / 2
    κ         = numerical_diffusivity(window_sum(budget.section_tendency), window_sum(budget.section_gradient))
    return (; label, budget.z, latitudes, κ)
end

plot_global_diffusivity_section(sections) #src
# ![Latitude–depth sections](global_diffusivity_section.png)
#
# WRK3-UP raises the diffusivity almost everywhere, while the differences of the other cases are small and change sign
# from one latitude to the next.

# ## Evolution in time
#
# Finally, we form the diffusivity of each 90-day window separately, for the whole ocean and for the levels below
# 1000 m, and place it at the middle of its window. This shows whether the differences between cases persist through
# the budget period or come from a few windows.

function diffusivity_in_layer(budget, levels)
    P = vec(sum(view(budget.wall_tendency, :, levels, :); dims=(1, 2)))
    G = vec(sum(view(budget.wall_gradient, :, levels, :); dims=(1, 2)))
    return numerical_diffusivity(P, G)
end

histories = map(labels, budgets) do label, budget
    return (; label, years=(budget.window_ends .- 45) / 365,
              κ=diffusivity_in_layer(budget, eachindex(budget.z)),
              deep_diffusivity=diffusivity_in_layer(budget, findall(≤(-1000), budget.z)))
end

plot_global_diffusivity_history(histories) #src
# ![Numerical diffusivity histories](global_diffusivity_history.png)
