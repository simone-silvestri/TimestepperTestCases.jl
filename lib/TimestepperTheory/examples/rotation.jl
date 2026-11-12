# # Coupled stability with rotation
# The boundaries compare frozen and progressive slow forcing over the outer scheme's wavenumber range.

include("rotation_scan.jl")
include("rotation_plots.jl")

function write_rotation_results(directory)
    results = run_scan(directory)
    plot_rotation_results(results, directory)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    write_rotation_results(joinpath(@__DIR__, "..", "output", "rotation"))
end
