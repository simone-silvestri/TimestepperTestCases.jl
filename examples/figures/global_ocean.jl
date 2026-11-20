# # Global ocean: computational cost
#
# The final timings of the quarter-degree configuration on one Nvidia H100 GPU are archived as `reported_timings.csv`
# in the Zenodo record. This example reads them, recomputes the simulated years per wall-clock day (SYPD) as
# Δt / (365 t_step) with both times in seconds, and the throughput relative to AB2-SE. Pass the directory holding
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

directory = isempty(ARGS) ? nothing : first(ARGS)
path      = dataset_file("global_ocean", "reported_timings.csv"; directory)
rows      = [split(line, ','; limit=5) for line in readlines(path)[2:end]]

measurements = map(rows) do row
    label, Δt, seconds_per_step, reported_years_per_day = row[1], parse(Float64, row[2]), parse(Float64, row[3]), parse(Float64, row[4])
    years_per_day = Δt / (365seconds_per_step)
    @assert round(years_per_day; digits=3) ≈ reported_years_per_day
    return (; label=String(label), Δt, seconds_per_step, years_per_day)
end

reference           = only(filter(measurement -> measurement.label == "AB2-SE", measurements))
relative_throughput = [measurement.years_per_day / reference.years_per_day for measurement in measurements]

for (measurement, throughput) in zip(measurements, relative_throughput)
    println(rpad(measurement.label, 14), round(measurement.years_per_day; digits=3), " SYPD, ", round(throughput; digits=2), " relative to AB2-SE")
end

plot_global_cost(measurements, relative_throughput) #src
# ![Global computational throughput](global_throughput.png)
