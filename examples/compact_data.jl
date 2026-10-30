using JLD2
using Statistics

function centered_output(file, name, iteration, grid)
    A = file["timeseries/$name/$iteration"]
    location = file["timeseries/$name/serialized/location"]
    N = (grid.Nx, grid.Ny, grid.Nz)
    H = (grid.Hx, grid.Hy, grid.Hz)
    for d in 1:3
        indices = ntuple(q -> q == d ? (H[d]+1:H[d]+N[d]) : Colon(), 3)
        if location[d] == Face && N[d] > 1
            neighbors = ntuple(q -> q == d ? (H[d]+2:H[d]+N[d]+1) : Colon(), 3)
            A = (A[indices...] + A[neighbors...]) / 2
        else
            A = A[indices...]
        end
    end
    return A
end

function prepare_channel(path)
    return jldopen(path) do file
        grid = file["serialized/grid"]
        iterations = sort(parse.(Int, keys(file["timeseries/t"])))
        iteration = last(iterations)
        fields = Dict(name => dropdims(mean(centered_output(file, name, iteration, grid); dims=1); dims=1)
                      for name in ("Abx", "Aby", "Abz", "Gbx", "Gby", "Gbz", "Dbz"))
        return (; y=collect(grid.yᵃᶜᵃ[1:grid.Ny]), z=collect(grid.z.cᵃᵃᶜ[1:grid.Nz]),
                  time=file["timeseries/t/$iteration"], fields)
    end
end

function prepare_internal_tide(path)
    return jldopen(path) do file
        saved_grid = file["serialized/grid"]
        grid = saved_grid.underlying_grid
        iterations = sort(parse.(Int, keys(file["timeseries/t"])))
        t = [file["timeseries/t/$iteration"] for iteration in iterations]
        z = collect(grid.z.cᵃᵃᶜ[1:grid.Nz])
        x = collect(grid.xᶜᵃᵃ[1:grid.Nx])
        bottom = vec(Array(saved_grid.immersed_boundary.bottom_height[1:grid.Nx, 1, 1]))
        wet = [depth > floor for floor in bottom, depth in z]
        fields = Dict(name => cat([dropdims(centered_output(file, name, iteration, grid); dims=2)
                                  for iteration in iterations]...; dims=3) for name in ("Abx", "Abz", "Gbz"))
        return (; x, z, Δz=diff(collect(grid.z.cᵃᵃᶠ[1:grid.Nz+1])), t, wet, fields)
    end
end

function compact_dataset(dataset, filename, prepare; directory=nothing)
    raw_directory = isnothing(directory) ? joinpath(@__DIR__, "..", dataset) : abspath(directory)
    raw_path = joinpath(raw_directory, filename)
    compact_name = "compact_" * filename
    local_compact = joinpath(raw_directory, compact_name)
    isfile(local_compact) && return JLD2.load(local_compact, "data")
    if isfile(raw_path)
        data = prepare(raw_path)
        output_directory = joinpath(@__DIR__, "..", "data", dataset)
        mkpath(output_directory)
        JLD2.jldsave(joinpath(output_directory, compact_name); compress=true, data)
        return data
    end
    return JLD2.load(dataset_file(dataset, compact_name; directory), "data")
end
