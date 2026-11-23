# # Global ocean: surface-average tracers
#
# In this example, we compare the surface temperature and salinity of the six global configurations over the five years, from
# the `global_ocean_<label>_average.jld2` files of the record. The figure is not part of the manuscript. Run it from
# the repository root with `julia --project=examples examples/figures/global_surface_averages.jl /path/to/global_ocean`.

using TimestepperTestCases
using Oceananigans
using JLD2
using Statistics
include("../shared/data_sources.jl") #src
include("../shared/compact_data.jl") #src
include("../shared/plotting.jl") #src
include("../shared/global_plots.jl") #src

directory = isempty(ARGS) ? joinpath(@__DIR__, "..", "..", "global_ocean") : first(ARGS)
labels = ["AB2-SE", "WRK3-SE-SM05", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-UP"]

# ## Area-weighted surface means
#
# The surface means check that every discretization keeps the same large-scale surface state while its eddies differ.
#
# `prepare_global_surface_averages` weights each surface cell by its area and keeps wet cells only. Times identify the
# ends of the saved averaging windows, and the record at time zero is the initial state. The means are cached in
# `data/global_ocean/`.

averages = map(labels) do label
    @info "Reading surface averages" label
    filename = "global_ocean_$(label)_average.jld2"
    compact_dataset("global_ocean", filename, prepare_global_surface_averages; directory, reuse=true)
end

# ## Differences from AB2-SE
#
# Every configuration saved its averages at the same times, so the differences from AB2-SE are taken record by record.

reference = first(averages)
@assert all(average -> average.days ≈ reference.days, averages)
histories = [(; label, years=average.days / 365, average.temperature, average.salinity,
                temperature_difference=average.temperature - reference.temperature,
                salinity_difference=average.salinity - reference.salinity) for (label, average) in zip(labels, averages)]

plot_global_surface_tracer_histories(histories) #src
# ![Area-averaged surface tracers](global_surface_tracer_histories.png)
#
# The Runge–Kutta runs stay within about 0.02 °C and 0.003 psu of AB2-SE over the five years. WRK3-UP drifts away: its
# surface cools by 0.1–0.15 °C and freshens slightly.
