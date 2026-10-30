# # Global configuration: computational cost
#
# This example reads the small timing records produced by global_ocean_cost. It uses the nominal production
# timestep to compute throughput, even when the timing experiment advances with a smaller cold-start timestep.
# Pass the directory containing global_ocean_cost_<label>_cost.jld2 as the first command-line argument.

using TimestepperTestCases
using JLD2
include("data_sources.jl") #src
include("plotting.jl") #src

directory = isempty(ARGS) ? nothing : first(ARGS)
measurements = NamedTuple[]

variants = global_ocean_discretizations()
if length(ARGS) > 1
    variants = [variant for variant in variants if variant.label in ARGS[2:end]]
end

for variant in variants
    filename = "global_ocean_cost_$(variant.label)_cost.jld2"
    path = dataset_file("global_ocean", filename; directory)
    record = JLD2.load(path)
    Δt = record["nominal_timestep"]
    seconds_per_step = record["seconds_per_step"]
    years_per_day = Δt / seconds_per_step / 365
    push!(measurements, (; label=variant.label, Δt, seconds_per_step, years_per_day))
end

reference = only(filter(measurement -> measurement.label == "AB2-SE", measurements))
relative_throughput = [measurement.years_per_day / reference.years_per_day for measurement in measurements]

for (measurement, throughput) in zip(measurements, relative_throughput)
    println(measurement.label, ": ", measurement.years_per_day, " simulated years/day; ", throughput, " relative to AB2-SE")
end

plot_global_cost(measurements, relative_throughput) #src
# ![Global computational throughput](global_throughput.png)
