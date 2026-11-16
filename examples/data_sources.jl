using Downloads: download
using SHA: sha256
using TOML: parsefile

function dataset_file(dataset, filename; directory=nothing)
    if !isnothing(directory)
        path = joinpath(abspath(directory), filename)
        isfile(path) || throw(ArgumentError("Local dataset file does not exist: $path"))
        return path
    end

    local_path = joinpath(@__DIR__, "..", dataset, filename)
    isfile(local_path) && return local_path

    record   = parsefile(joinpath(@__DIR__, "datasets.toml"))[dataset]
    checksum = get(record["sha256"], filename, "")
    path     = joinpath(@__DIR__, "..", "data", dataset, filename)
    if isfile(path) && !isempty(checksum) && open(sha256, path) == hex2bytes(checksum)
        return path
    end
    isempty(record["base_url"]) && throw(ArgumentError("Dataset $dataset is not published yet; supply a local data directory."))
    isempty(checksum) && throw(ArgumentError("The download checksum for $filename is missing from datasets.toml."))
    mkpath(dirname(path))
    mktemp(dirname(path)) do temporary_path, stream
        close(stream)
        download(rstrip(record["base_url"], '/') * "/" * filename, temporary_path)
        open(sha256, temporary_path) == hex2bytes(checksum) || error("Checksum mismatch for $filename.")
        mv(temporary_path, path; force=true)
    end
    return path
end
