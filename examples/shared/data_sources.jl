include("../../scripts/download_data.jl")
using TOML: parsefile

function dataset_file(dataset, filename; directory=nothing)
    if !isnothing(directory)
        path = joinpath(abspath(directory), filename)
        isfile(path) || throw(ArgumentError("Local dataset file does not exist: $path"))
        return path
    end

    local_path = joinpath(@__DIR__, "..", "..", dataset, filename)
    isfile(local_path) && return local_path

    record   = parsefile(joinpath(@__DIR__, "..", "datasets.toml"))[dataset]
    checksum = get(record["sha256"], filename, "")
    path     = joinpath(@__DIR__, "..", "..", "data", dataset, filename)
    if isfile(path) && !isempty(checksum) && open(sha256, path) == hex2bytes(checksum)
        return path
    end
    isempty(record["base_url"]) && throw(ArgumentError("Dataset $dataset is not published yet; supply a local data directory."))
    isempty(checksum) && throw(ArgumentError("The download checksum for $filename is missing from datasets.toml."))
    return download_verified(rstrip(record["base_url"], '/') * "/" * filename, path, checksum)
end
