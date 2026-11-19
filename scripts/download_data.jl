using Downloads: download
using SHA: sha256

function download_verified(url, path, checksum)
    expected = hex2bytes(checksum)
    isfile(path) && open(sha256, path) == expected && return path
    mkpath(dirname(abspath(path)))
    mktemp(dirname(abspath(path))) do temporary_path, stream
        close(stream)
        download(url, temporary_path)
        open(sha256, temporary_path) == expected || error("Checksum mismatch for $(basename(path)).")
        mv(temporary_path, path; force=true)
    end
    return path
end

function download_campaigns(record, directory)
    record > 0 || throw(ArgumentError("Supply a positive, version-specific Zenodo record ID."))
    directory = abspath(directory)
    base_url = "https://zenodo.org/records/$record/files/"
    mkpath(directory)
    checksums = download(base_url * "SHA256SUMS")
    try
        for line in eachline(checksums)
            isempty(strip(line)) && continue
            checksum, relative_path = split(line; limit=2)
            relative_path = strip(relative_path)
            components = splitpath(relative_path)
            (isabspath(relative_path) || ".." in components) && error("Invalid archive path: $relative_path")
            path = joinpath(directory, relative_path)
            @info "Checking dataset file" path
            download_verified(base_url * basename(path), path, checksum)
        end
        cp(checksums, joinpath(directory, "SHA256SUMS"); force=true)
    finally
        rm(checksums; force=true)
    end
    return directory
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2 || first(ARGS) in ("-h", "--help")
        println("Usage: julia scripts/download_data.jl RECORD_ID DESTINATION")
        println("Restore the three campaign folders and verify every file against SHA256SUMS.")
    else
        download_campaigns(parse(Int, ARGS[1]), ARGS[2])
    end
end
