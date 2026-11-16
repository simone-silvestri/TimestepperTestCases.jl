# # Global ocean: numerical diffusivity
#
# Run `julia --project=examples examples/global_ocean_results.jl global_ocean/raw` to reproduce these figures.
# The six campaign files are named `global_ocean_<label>_dissipation.jld2`. Compact reductions are saved in
# `data/global_ocean/`; pass that directory to redraw the figures without the original three-dimensional fields.

using TimestepperTestCases
using Oceananigans
using Statistics
include("data_sources.jl") #src
include("compact_data.jl") #src
include("plotting.jl") #src

directory = isempty(ARGS) ? joinpath(@__DIR__, "..", "global_ocean", "raw") : first(ARGS)
labels = ["AB2-SE", "WRK3-SE-SM05", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-UP"]
budgets = map(labels) do label
    @info "Reading global buoyancy-variance budget" label
    filename = "global_ocean_$(label)_dissipation.jld2"
    compact_dataset("global_ocean", filename, prepare_global_budget; directory, reuse=true)
end

# ## Numerical diffusivity
#
# Form κ = −ΣP/(2ΣG) after integrating the variance tendency P and squared buoyancy gradient G over each region.
# Both quantities already include cell volumes. The reduction places staggered diagnostics at cell centers and
# excludes nonfinite metric values from both integrals. The first budget output straddles the diagnostic startup;
# the eight complete 90-day windows end at days 1170 through 1800. Spatial averages use these same windows in all cases.

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

plot_global_diffusivity_wall_distance(profiles) #src
# ![Profiles by distance from topography](global_diffusivity_wall_distance.png)

# The first contour panel shows WRK3-SE on a logarithmic scale from 10⁻⁶ to 10⁻³ m² s⁻¹; values below this range
# use the lowest color. The remaining panels show signed differences from WRK3-SE on a separate diverging scale.

sections = map(labels, budgets) do label, budget
    latitudes = (budget.latitude_edges[1:end-1] + budget.latitude_edges[2:end]) / 2
    κ         = numerical_diffusivity(window_sum(budget.section_tendency), window_sum(budget.section_gradient))
    return (; label, budget.z, latitudes, κ)
end

plot_global_diffusivity_section(sections) #src
# ![Latitude–depth sections](global_diffusivity_section.png)

# ## Evolution over the complete budget windows
#
# Each point is placed at the midpoint of its 90-day window. The deep-ocean integral uses levels below 1000 m.

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
