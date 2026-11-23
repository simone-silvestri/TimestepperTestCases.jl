# # Global ocean: computational cost
#
# The final timings of the quarter-degree configuration on one Nvidia H100 GPU are archived as `reported_timings.csv`
# in the Zenodo record. In this example, we read them, recompute the simulated years per wall-clock day (SYPD)
# as Δt / (365 t_step) with both times in seconds, and compute the throughput relative to AB2-SE. Pass the directory holding
# `reported_timings.csv` as the first argument; without it, the file is downloaded from the record and verified.
#
# | Configuration | Δt (min) | t_step (s) | SYPD   | Relative throughput |
# |:--------------|---------:|-----------:|-------:|--------------------:|
# | AB2-SE        | 11       | 0.2534     | 7.136  | 1.00                |
# | WRK3-SE       | 45       | 0.9097     | 8.132  | 1.14                |
# | SRK3-SE       | 45       | 0.8873     | 8.337  | 1.17                |
# | MRK4-SE       | 72       | 1.3859     | 8.540  | 1.20                |
# | WRK3-SE-F1    | 45       | 0.8284     | 8.930  | 1.25                |
# | WRK3-IM       | 45       | 2.1238     | 3.483  | 0.49                |
# | WRK3-UP       | 45       | 0.7386     | 10.015 | 1.40                |
#
# WRK3-SE-F1 is labelled WRK3-SE-SM05 in the code and file names. WRK3-IM was timed but not integrated for five years.

using TimestepperTestCases
include("../shared/data_sources.jl") #src
include("../shared/plotting.jl") #src
include("../shared/global_plots.jl") #src

# ## Reading the timings
#
# Each row of `reported_timings.csv` holds a configuration, its time step and its measured wall time per step, both in
# seconds, together with the reported SYPD.

directory = isempty(ARGS) ? nothing : first(ARGS)
path      = dataset_file("global_ocean", "reported_timings.csv"; directory)
rows      = [split(line, ','; limit=5) for line in readlines(path)[2:end]]

# ## Throughput
#
# A time step Δt that costs t_step seconds of wall time advances the model by Δt / t_step simulated days per wall-clock
# day, hence Δt / (365 t_step) simulated years per day. The recomputed value must match the reported one to the three
# digits of the table.

measurements = map(rows) do row
    label, Δt, seconds_per_step, reported_years_per_day = row[1], parse(Float64, row[2]), parse(Float64, row[3]), parse(Float64, row[4])
    years_per_day = Δt / (365seconds_per_step)
    @assert round(years_per_day; digits=3) ≈ reported_years_per_day
    return (; label=String(label), Δt, seconds_per_step, years_per_day)
end

# The throughput relative to AB2-SE is the ratio of the SYPD.

reference           = only(filter(measurement -> measurement.label == "AB2-SE", measurements))
relative_throughput = [measurement.years_per_day / reference.years_per_day for measurement in measurements]

for (measurement, throughput) in zip(measurements, relative_throughput)
    println(rpad(measurement.label, 14), round(measurement.years_per_day; digits=3), " SYPD, ", round(throughput; digits=2), " relative to AB2-SE")
end

plot_global_cost(measurements, relative_throughput) #src
# ![Global computational throughput](global_throughput.png)
#
# A three-stage Runge–Kutta step costs about 3.5 times an AB2 step, but it is about four times longer, so WRK3-SE,
# SRK3-SE and MRK4-SE run 14 to 20% faster than AB2-SE. The barotropic treatment matters less: WRK3-SE-F1, with about
# a quarter of the barotropic work of WRK3-SE, is only 10% faster. The implicit free surface halves the throughput,
# because its global elliptic solve is expensive on a GPU.
