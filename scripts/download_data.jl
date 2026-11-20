using Downloads: download, request, RequestError
using SHA: sha256

const zenodo_record       = 24315107
const published_campaigns = ("internal_tide", "channel", "global_ocean")

"""
    download_verified(url, path, checksum; attempts=10)

Download `url` to `path` and verify its SHA-256 `checksum`. The transfer accumulates in `path * ".partial"`, so
an interrupted download resumes from the received bytes, both within the `attempts` and when rerun.
"""
function download_verified(url, path, checksum; attempts=10)
    expected = hex2bytes(checksum)
    isfile(path) && open(sha256, path) == expected && return path
    mkpath(dirname(abspath(path)))
    partial_path = path * ".partial"
    for attempt in 1:attempts
        offset   = filesize(partial_path)
        headers  = offset > 0 ? ["Range" => "bytes=$offset-"] : Pair{String, String}[]
        response = try
            open(stream -> request(url; output=stream, headers), partial_path, "a")
        catch exception
            exception isa RequestError || rethrow()
            exception
        end
        status = response isa RequestError ? 0 : response.status
        status in (0, 200, 206) || open(stream -> truncate(stream, offset), partial_path, "r+") # drop the error body
        if status == 404
            rm(partial_path; force=true)
            throw(RequestError(url, 0, "Not found", response))
        elseif status == 200 && offset > 0
            rm(partial_path) # the server ignored the range and appended the whole file
        elseif status in (200, 206, 416)
            break
        end
        attempt == attempts && error("Could not download $(basename(path)) in $attempts attempts; rerun to resume.")
        @warn "Resuming interrupted download" file=basename(path) received_bytes=filesize(partial_path) attempt
        sleep(10)
    end
    if open(sha256, partial_path) != expected
        rm(partial_path)
        error("Checksum mismatch for $(basename(path)); rerun to download it again.")
    end
    chmod(partial_path, 0o644)
    mv(partial_path, path; force=true)
    return path
end

function archive_entries(checksums)
    lines = filter(line -> !isempty(strip(line)), readlines(checksums))
    return map(lines) do line
        checksum, relative_path = split(line; limit=2)
        relative_path = strip(relative_path)
        components = splitpath(relative_path)
        (isabspath(relative_path) || ".." in components) && error("Invalid archive path: $relative_path")
        return (; checksum, relative_path, campaign=length(components) == 1 ? nothing : first(components))
    end
end

"""
    download_campaigns(directory, campaigns=published_campaigns; record=zenodo_record, concurrent_downloads=4)

Download the files of `campaigns` from the Zenodo `record` into `directory`, restore the campaign folders, and
verify every file against the record's `SHA256SUMS`. Files already present with a matching checksum are reused.
Checksum entries without a file in the record are reported at the end.
"""
function download_campaigns(directory, campaigns=published_campaigns; record=zenodo_record, concurrent_downloads=4)
    record > 0 || throw(ArgumentError("Supply a positive, version-specific Zenodo record ID."))
    unknown = setdiff(campaigns, published_campaigns)
    isempty(unknown) || throw(ArgumentError("Unknown campaigns $(join(unknown, ", ")); choose from $(join(published_campaigns, ", "))."))
    directory = abspath(directory)
    base_url = "https://zenodo.org/records/$record/files/"
    mkpath(directory)
    checksums = download(base_url * "SHA256SUMS")
    try
        entries = filter(entry -> isnothing(entry.campaign) || entry.campaign in campaigns, archive_entries(checksums))
        absent_files = asyncmap(entries; ntasks=concurrent_downloads) do entry
            path = joinpath(directory, entry.relative_path)
            @info "Checking dataset file" path
            try
                download_verified(base_url * basename(path), path, entry.checksum)
                return nothing
            catch exception
                exception isa RequestError && exception.response.status == 404 || rethrow()
                return entry.relative_path
            end
        end
        filter!(!isnothing, absent_files)
        isempty(absent_files) || @info "Skipped files that SHA256SUMS lists but record $record does not provide" absent_files
        cp(checksums, joinpath(directory, "SHA256SUMS"); force=true)
    finally
        rm(checksums; force=true)
    end
    return directory
end

if abspath(PROGRAM_FILE) == @__FILE__
    arguments = copy(ARGS)
    record    = !isempty(arguments) && all(isdigit, first(arguments)) ? parse(Int, popfirst!(arguments)) : zenodo_record
    if isempty(arguments) || first(arguments) in ("-h", "--help")
        println("Usage: julia scripts/download_data.jl [RECORD_ID] DESTINATION [CAMPAIGN ...]")
        println("Download the campaigns ($(join(published_campaigns, ", ")); all by default) from Zenodo record")
        println("$zenodo_record, restore their folders, and verify every file against SHA256SUMS.")
    else
        destination = popfirst!(arguments)
        download_campaigns(destination, isempty(arguments) ? published_campaigns : Tuple(arguments); record)
    end
end
