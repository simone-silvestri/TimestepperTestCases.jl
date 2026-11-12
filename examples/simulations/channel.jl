using CUDA
using TimestepperTestCases
using Oceananigans

# Usage: julia --project=. channel.jl [labels...]
#
# With no labels every case runs in turn. If `CHANNEL_RESTART_DIR` is set, each case starts from its spin-up
# restart in that directory instead of repeating the spin-up. Output is written to the working directory.

# Redirected output is buffered, so it is flushed every minute to keep a log current
flusher = Timer(_ -> (flush(stdout); flush(stderr)), 60; interval = 60)

labels            = isempty(ARGS) ? [d.label for d in discretizations()] : ARGS
restart_directory = get(ENV, "CHANNEL_RESTART_DIR", nothing)

for label in labels
    restart_file = isnothing(restart_directory) ? nothing : joinpath(restart_directory, "restart$(label)_iteration144000.jld2")
    channel_simulation(discretization(label); restart_file)
end
