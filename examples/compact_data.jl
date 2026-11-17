using JLD2
using Statistics
using TimestepperTestCases: power_spectrum_x

function centered_output(file, name, iteration, grid)
    A        = file["timeseries/$name/$iteration"]
    location = file["timeseries/$name/serialized/location"]
    N        = (grid.Nx, grid.Ny, grid.Nz)
    H        = (grid.Hx, grid.Hy, grid.Hz)
    for d in 1:3
        indices = ntuple(q -> q == d ? (H[d]+1:H[d]+N[d]) : Colon(), 3)
        if location[d] == Face && N[d] > 1
            neighbors = ntuple(q -> q == d ? (H[d]+2:H[d]+N[d]+1) : Colon(), 3)
            A         = (A[indices...] + A[neighbors...]) / 2
        else
            A = A[indices...]
        end
    end
    return A
end

function prepare_channel(path)
    return jldopen(path) do file
        grid       = file["serialized/grid"]
        iterations = sort(parse.(Int, keys(file["timeseries/t"])))
        iteration  = last(iterations)
        fields     = Dict(name => dropdims(mean(centered_output(file, name, iteration, grid); dims=1); dims=1)
                      for name in ("Abx", "Aby", "Abz", "Gbx", "Gby", "Gbz", "Dbz"))
        return (; y=collect(grid.yᵃᶜᵃ[1:grid.Ny]), z=collect(grid.z.cᵃᵃᶜ[1:grid.Nz]),
                  time=file["timeseries/t/$iteration"], fields)
    end
end

function prepare_internal_tide(path)
    return jldopen(path) do file
        saved_grid = file["serialized/grid"]
        grid       = saved_grid.underlying_grid
        iterations = sort(parse.(Int, keys(file["timeseries/t"])))
        t          = [file["timeseries/t/$iteration"] for iteration in iterations]
        z          = collect(grid.z.cᵃᵃᶜ[1:grid.Nz])
        x          = collect(grid.xᶜᵃᵃ[1:grid.Nx])
        bottom     = vec(Array(saved_grid.immersed_boundary.bottom_height[1:grid.Nx, 1, 1]))
        wet        = [depth > floor for floor in bottom, depth in z]
        fields     = Dict(name => cat([dropdims(centered_output(file, name, iteration, grid); dims=2)
                                  for iteration in iterations]...; dims=3) for name in ("Abx", "Abz", "Gbz"))
        return (; x, z, Δz=diff(collect(grid.z.cᵃᵃᶠ[1:grid.Nz+1])), t, wet, fields)
    end
end

#####
##### Global configuration
#####

function global_ocean_geometry(path)
    saved_grid = jldopen(file -> file["serialized/grid"], path)
    grid       = saved_grid.underlying_grid
    (; Nx, Ny, Nz, Hx, Hy, Hz) = grid
    centers(name) = Array(getproperty(grid, name)[1:Nx, 1:Ny])
    corners(name) = Array(getproperty(grid, name)[1:Nx+1, 1:Ny+1])
    return (; Nx, Ny, Nz, Hx, Hy, Hz, wet=wet_mask(saved_grid.interior_active_cells, Nx, Ny, Nz),
              λ=centers(:λᶜᶜᵃ), φ=centers(:φᶜᶜᵃ), Δx=centers(:Δxᶜᶜᵃ), Az=centers(:Azᶜᶜᵃ),
              λᶠᶠ=corners(:λᶠᶠᵃ), φᶠᶠ=corners(:φᶠᶠᵃ),
              zᶜ=collect(grid.z.cᵃᵃᶜ[1:Nz]), zᶠ=collect(grid.z.cᵃᵃᶠ[1:Nz+1]))
end

function wet_mask(active_cells, N...)
    wet = falses(N...)
    for cell in active_cells
        wet[cell...] = true
    end
    return wet
end

latitude_bands(φ, latitude_edges) = [clamp(searchsortedlast(latitude_edges, value), 1, length(latitude_edges) - 1) for value in φ]

# Restarts leave extra records shortly after a whole day, so each day keeps only its closest record.
function daily_records(file)
    iterations = sort(parse.(Int, keys(file["timeseries/t"])))
    days       = [file["timeseries/t/$iteration"] / 86400 for iteration in iterations]
    closest    = [argmin(abs.(days .- day)) for day in round(Int, first(days)):round(Int, last(days))]
    return iterations[closest], days[closest]
end

function surface_kinetic_energy(u, v, geometry)
    (; Nx, Ny, Nz, Hx, Hy, wet) = geometry
    return [wet[i, j, Nz] ? (u[Hx+i, Hy+j, 1]^2 + u[Hx+mod1(i+1, Nx), Hy+j, 1]^2 + v[Hx+i, Hy+j, 1]^2 + v[Hx+i, Hy+j+1, 1]^2) / 4 : NaN32
            for i in 1:Nx, j in 1:Ny]
end

# Stretches of `segment_length` cells along the rows of a region, overlapping by three quarters. A cell enters only if the
# velocity faces on its western and southern sides separate two wet cells.
function land_free_segments(region, geometry; segment_length)
    (; Nx, Ny, Nz, wet, λ, φ) = geometry
    inside(i, j) = first(region.longitudes) ≤ mod(λ[i, j], 360) ≤ last(region.longitudes) && first(region.latitudes) ≤ φ[i, j] ≤ last(region.latitudes)
    open_water(i, j) = j > 1 && inside(i, j) && wet[i, j, Nz] && wet[mod1(i - 1, Nx), j, Nz] && wet[i, j-1, Nz]
    columns(i) = mod1.(i:i+segment_length-1, Nx)
    return [(; row=j, columns=columns(i)) for j in 1:Ny for i in 1:segment_length÷4:Nx if all(column -> open_water(column, j), columns(i))]
end

# Spectral density of a non-periodic segment: the linear trend is removed and the Hann window has a unit mean square, so that
# the wavenumber integral returns the variance of the detrended signal.
function tapered_spectrum(values, Δx)
    L         = length(values)
    n         = eachindex(values) .- (L + 1) / 2
    detrended = values .- mean(values) .- n * sum(n .* values) / sum(abs2, n)
    taper     = @. sqrt(8 / 3) * cos(π * n / L)^2
    return real.(first(power_spectrum_x(taper .* detrended, (0:L-1) .* Δx))) .* (L * Δx / 2π)
end

# Cells to the nearest land along a grid line, the only direction a WENO stencil reaches. The rigid bottom and the
# southern edge count as land; the free surface and the northern fold do not.
function wall_distance(wet; reach=5)
    Nx, Ny, Nz = size(wet)
    wall_at(i, j, k, n) = !wet[mod1(i + n, Nx), j, k] || !wet[mod1(i - n, Nx), j, k] ||
                          (j + n ≤ Ny && !wet[i, j+n, k]) || j ≤ n || !wet[i, j-n, k] ||
                          (k + n ≤ Nz && !wet[i, j, k+n]) || k ≤ n || !wet[i, j, k-n]
    distance = fill(Int8(reach + 1), size(wet))
    Threads.@threads for k in 1:Nz
        for j in 1:Ny, i in 1:Nx
            wet[i, j, k] || continue
            n = findfirst(n -> wall_at(i, j, k, n), 1:reach)
            isnothing(n) || (distance[i, j, k] = n)
        end
    end
    return distance
end

# Dry cells are skipped because the gradients carry NaN on some of their faces, and the faces are summed in Float64.
function centered_total(x, y, z, geometry)
    (; Nx, Ny, Nz, Hx, Hy, Hz, wet) = geometry
    total = zeros(Nx, Ny, Nz)
    Threads.@threads for k in 1:Nz
        for j in 1:Ny, i in 1:Nx
            wet[i, j, k] || continue
            I, J, K = Hx + i, Hy + j, Hz + k
            total[i, j, k] = (Float64(x[I, J, K]) + x[Hx+mod1(i+1, Nx), J, K] + y[I, J, K] + y[I, J+1, K] + z[I, J, K] + z[I, J, K+1]) / 2
        end
    end
    return total
end

read_centered_total(file, names, iteration, geometry) = centered_total((file["timeseries/$name/$iteration"] for name in names)..., geometry)

function binned_level_sum(values, wet, bin, number_of_bins)
    Nx, Ny, Nz = size(values)
    totals = zeros(number_of_bins, Nz)
    Threads.@threads for k in 1:Nz
        for j in 1:Ny, i in 1:Nx
            wet[i, j, k] && (totals[bin(i, j, k), k] += values[i, j, k])
        end
    end
    return totals
end

function prepare_global_budget(path; latitude_edges=-80:0.25:90, minimum_day=1170)
    geometry = global_ocean_geometry(path)
    (; wet, φ, zᶜ) = geometry
    distance            = wall_distance(wet)
    bands               = latitude_bands(φ, latitude_edges)
    number_of_distances = Int(maximum(distance))
    number_of_bands     = length(latitude_edges) - 1
    distance_bin(i, j, k) = distance[i, j, k]
    band_bin(i, j, k) = bands[i, j]

    return jldopen(path) do file
        iterations = sort(parse.(Int, keys(file["timeseries/t"])))
        filter!(iteration -> file["timeseries/t/$iteration"] / 86400 ≥ minimum_day, iterations)
        window_ends = [file["timeseries/t/$iteration"] / 86400 for iteration in iterations]
        windows     = map(iterations) do iteration
            P       = read_centered_total(file, ("Abx", "Aby", "Abz"), iteration, geometry)
            G       = read_centered_total(file, ("Gbx", "Gby", "Gbz"), iteration, geometry)
            invalid = .!isfinite.(P) .| .!isfinite.(G)
            P[invalid] .= 0
            G[invalid] .= 0
            @info "Reduced budget window" day=file["timeseries/t/$iteration"] / 86400 excluded_cells=count(invalid)
            (; wall_tendency=binned_level_sum(P, wet, distance_bin, number_of_distances),
               wall_gradient=binned_level_sum(G, wet, distance_bin, number_of_distances),
               section_tendency=binned_level_sum(P, wet, band_bin, number_of_bands),
               section_gradient=binned_level_sum(G, wet, band_bin, number_of_bands))
        end
        stacked(name) = stack(getproperty(window, name) for window in windows)
        return (; z=zᶜ, latitude_edges=collect(latitude_edges), window_ends,
                  wall_tendency=stacked(:wall_tendency), wall_gradient=stacked(:wall_gradient),
                  section_tendency=stacked(:section_tendency), section_gradient=stacked(:section_gradient))
    end
end

#####
##### Compact files
#####

# With `reuse`, a compact file already in the data cache spares the pass over the original output.
function compact_dataset(dataset, filename, prepare; directory=nothing, reuse=false)
    raw_directory = isnothing(directory) ? joinpath(@__DIR__, "..", dataset) : abspath(directory)
    raw_path      = joinpath(raw_directory, filename)
    compact_name  = "compact_" * filename
    local_compact = joinpath(raw_directory, compact_name)
    isfile(local_compact) && return JLD2.load(local_compact, "data")
    cached_compact = joinpath(@__DIR__, "..", "data", dataset, compact_name)
    if reuse && isfile(cached_compact) && (!isfile(raw_path) || mtime(cached_compact) ≥ mtime(raw_path))
        return JLD2.load(cached_compact, "data")
    end
    if isfile(raw_path)
        data = prepare(raw_path)
        mkpath(dirname(cached_compact))
        JLD2.jldsave(cached_compact; compress=true, data)
        return data
    end
    return JLD2.load(dataset_file(dataset, compact_name; directory), "data")
end

function prepare_global_surface_averages(path)
    (; Nx, Ny, Nz, Hx, Hy, wet, Az) = global_ocean_geometry(path)
    weights = Az .* view(wet, :, :, Nz)
    area    = sum(weights)
    return jldopen(path) do file
        iterations = sort(parse.(Int, keys(file["timeseries/t"])))
        days       = [file["timeseries/t/$iteration"] / 86400 for iteration in iterations]
        surface_mean(name, iteration) = sum(view(file["timeseries/$name/$iteration"], Hx+1:Hx+Nx, Hy+1:Hy+Ny, 1) .* weights) / area
        return (; days, temperature=[surface_mean("T", iteration) for iteration in iterations],
                  salinity=[surface_mean("S", iteration) for iteration in iterations])
    end
end
