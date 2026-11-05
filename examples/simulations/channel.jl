using CUDA
using TimestepperTestCases
using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces

# The cases of Table 1, plus the SM05 reference and the four-stage composition of section 4, defined once in
# `discretizations()` and shared by every test case.
#
#   WRK3-SE, WRK3-SE-SM05, WRK3-IM, SRK3-SE, AB2-SE, WRK3-UP, MRK4-SE
#
# The channel's stratification is surface-intensified, so its c₁ comes from the WKB integral -- see
# `channel_stability_parameters`.
#
# Usage: julia --project=. channel.jl [labels...]
#
# With no labels every case runs in turn. If `CHANNEL_RESTART_DIR` is set, each case starts from its spin-up
# restart in that directory instead of repeating the spin-up. Output is written to the working directory.

# Redirected output is buffered, so it is flushed every minute to keep a log current
flusher = Timer(_ -> (flush(stdout); flush(stderr)), 60; interval = 60)

labels      = isempty(ARGS) ? [d.label for d in discretizations()] : ARGS
restart_dir = get(ENV, "CHANNEL_RESTART_DIR", nothing)

for label in labels
    restart_file = isnothing(restart_dir) ? nothing : joinpath(restart_dir, "restart$(label)_iteration144000.jld2")
    sim = channel_simulation(discretization(label); restart_file)
end
