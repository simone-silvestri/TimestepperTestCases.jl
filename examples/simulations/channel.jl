# # Channel campaign
#
# In this script, we run the re-entrant channel for every discretization, or only for the labels passed as arguments:
# `julia --project=. examples/simulations/channel.jl [labels...]`. Each run spins up for 100 days at a one-minute step
# and then integrates for 40 years; the output is written to the working directory.

using CUDA
using TimestepperTestCases
using Oceananigans

# With `CHANNEL_RESTART_DIR` set, each case starts from its spin-up restart in that directory instead of repeating the
# spin-up.

# Redirected output is buffered, so it is flushed every minute to keep a log current
flusher = Timer(_ -> (flush(stdout); flush(stderr)), 60; interval = 60)

labels            = isempty(ARGS) ? [d.label for d in discretizations()] : ARGS
restart_directory = get(ENV, "CHANNEL_RESTART_DIR", nothing)

for label in labels
    restart_file = isnothing(restart_directory) ? nothing : joinpath(restart_directory, "restart$(label)_iteration144000.jld2")
    channel_simulation(discretization(label); restart_file)
end
