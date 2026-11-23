# # AB2-SE and the output schedule
#
# In this script, we run AB2-SE twice, saving output every 18 steps and every hour. Output saved at fixed times shortens
# the step that lands on each output time, and the quasi-AB2 scheme restarts a shortened step with forward Euler. We
# record the time step, the restarts and the horizontally averaged variance budget at every step of both runs.

using TimestepperTestCases, Oceananigans, JLD2, Statistics
using Oceananigans.Units
using Oceananigans.Fields: location

root = abspath(get(ENV, "AB2_SCHEDULE_OUTPUT", joinpath(@__DIR__, "..", "..", "internal_tide", "schedule_experiment")))
mkpath(root)

function centered_profile(field)
    A = interior(field)
    if location(field, 1) == Face && size(A, 1) > 1
        A = (A + circshift(A, (-1, 0, 0))) / 2
    end
    if location(field, 3) == Face && size(A, 3) > 1
        A = (A[:, :, 1:end-1] + A[:, :, 2:end]) / 2
    end
    return vec(sum(A; dims=(1, 2)))
end

for experiment in ("iteration", "time")
    directory = joinpath(root, experiment)
    mkpath(directory)
    cd(directory) do
        schedule = experiment == "iteration" ? IterationInterval(18) : TimeInterval(1hours)
        simulation = internal_tide(discretization("AB2-SE"); output_schedule=schedule, run_simulation=false)
        outputs = simulation.output_writers[:fields].outputs
        gradient = Field(outputs.Gbz)
        times = Float64[]
        steps = Float64[]
        restarts = Bool[]
        horizontal = Vector{Float64}[]
        vertical = Vector{Float64}[]
        gradients = Vector{Float64}[]
        previous_step = Ref(Inf)
        function record_step(sim)
            sim.model.clock.iteration == 0 && return nothing
            step = sim.model.clock.last_Δt
            compute!(gradient)
            push!(times, sim.model.clock.time)
            push!(steps, step)
            push!(restarts, step != previous_step[])
            push!(horizontal, centered_profile(outputs.Abx))
            push!(vertical, centered_profile(outputs.Abz))
            push!(gradients, centered_profile(gradient))
            previous_step[] = step
            return nothing
        end
        simulation.callbacks[:schedule_diagnostics] = Callback(record_step, IterationInterval(1))
        @info "Schedule experiment" experiment timestep=simulation.Δt schedule
        run!(simulation)
        jldsave("every_step.jld2"; times, steps, restarts, horizontal=stack(horizontal), vertical=stack(vertical),
                gradients=stack(gradients), nominal_timestep=simulation.Δt,
                oceananigans_source=pathof(Oceananigans))
        @info "Finished schedule experiment" experiment steps=length(times) restarts=count(restarts)
    end
    GC.gc()
end
