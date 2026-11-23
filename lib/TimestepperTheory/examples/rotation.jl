# # Coupled stability with rotation
#
# In this example, we extend the coupled stability analysis to a rotating system (Figure 19 of the manuscript). The
# stability boundaries compare frozen and progressive slow forcing over the wavenumber range of each outer scheme. Run it
# with `julia --project=lib/TimestepperTheory/examples --threads=auto lib/TimestepperTheory/examples/rotation.jl`.

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
