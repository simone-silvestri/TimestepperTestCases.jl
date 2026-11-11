# # Global ocean: surface-average tracers
#
# Run `julia --project=examples examples/global_surface_averages.jl global_ocean/global` from the repository root.
# These histories use the surface temperature and salinity in the six `_average.jld2` files. Each spatial mean
# is weighted by wet surface-cell area. Times identify the ends of the saved averaging windows; the record at
# zero is the initial state. The lower panels show differences from AB2-SE at the same output times.

using TimestepperTestCases
using Oceananigans
using JLD2
using Statistics
include("data_sources.jl") #src
include("compact_data.jl") #src
include("plotting.jl") #src

directory = isempty(ARGS) ? joinpath(@__DIR__, "..", "global_ocean", "global") : first(ARGS)
labels = ["AB2-SE", "WRK3-SE-SM05", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-UP"]
averages = map(labels) do label
    @info "Reading surface averages" label
    filename = "global_ocean_$(label)_average.jld2"
    compact_dataset("global_ocean", filename, prepare_global_surface_averages; directory, reuse=true)
end

reference = first(averages)
@assert all(average -> average.days ≈ reference.days, averages)
histories = [(; label, years=average.days / 365, average.temperature, average.salinity,
                temperature_difference=average.temperature - reference.temperature,
                salinity_difference=average.salinity - reference.salinity) for (label, average) in zip(labels, averages)]

plot_global_surface_tracer_histories(histories) #src
# ![Area-averaged surface tracers](global_surface_tracer_histories.png)
