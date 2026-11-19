# Retain daily surface velocities over the final year, including Oceananigans grid and field metadata.
# Run `julia --project=examples scripts/archive_global_surface.jl INPUT_DIRECTORY OUTPUT_DIRECTORY`.

using JLD2
using Oceananigans
include("../examples/shared/compact_data.jl")

function copy_output_group!(destination, source)
    for name in keys(source)
        value = source[name]
        if value isa JLD2.Group
            copy_output_group!(JLD2.Group(destination, name), value)
        else
            destination[name] = value
        end
    end
    return nothing
end

function archive_global_surface(source, destination)
    abspath(source) == abspath(destination) && error("Source and destination must differ.")
    ispath(destination) && error("Destination exists: $destination")
    temporary = destination * ".partial.jld2"
    ispath(temporary) && error("Incomplete archive exists: $temporary")
    mkpath(dirname(destination))
    jldopen(source) do input
        iterations, days = daily_records(input)
        retained = unique(iterations[findall(≥(last(days) - 365), days)])
        jldopen(temporary, "w"; compress=[JLD2.Shuffle(), JLD2.ZstdFilter(1)]) do output
            for name in keys(input)
                name == "timeseries" && continue
                value = input[name]
                value isa JLD2.Group ? copy_output_group!(JLD2.Group(output, name), value) : (output[name] = value)
            end
            for name in ("u", "v")
                copy_output_group!(JLD2.Group(output, "timeseries/$name/serialized"), input["timeseries/$name/serialized"])
                for iteration in retained
                    key = "timeseries/$name/$iteration"
                    output[key] = input[key]
                end
            end
            for iteration in retained
                key = "timeseries/t/$iteration"
                output[key] = input[key]
            end
        end
        jldopen(temporary) do output
            archived_iterations, archived_days = daily_records(output)
            selected = findall(≥(last(days) - 365), days)[1:5:end]
            @assert archived_iterations[1:5:end] == iterations[selected]
            @assert archived_days[1:5:end] == days[selected]
            @assert last(archived_iterations) == last(iterations)
            for name in ("u", "v"), iteration in retained
                key = "timeseries/$name/$iteration"
                @assert isequal(output[key], input[key])
            end
        end
    end
    for name in ("u", "v")
        series = FieldTimeSeries(temporary, name; backend=OnDisk())
        @assert !isempty(interior(series[1]))
        @assert !isempty(interior(series[length(series.times)]))
    end
    mv(temporary, destination)
    @info "Archived surface velocities" source destination bytes=filesize(destination)
    return destination
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error("Provide input and output directories.")
    for filename in sort(readdir(ARGS[1]))
        endswith(filename, "_surface.jld2") || continue
        archive_global_surface(joinpath(ARGS[1], filename), joinpath(ARGS[2], filename))
    end
end
