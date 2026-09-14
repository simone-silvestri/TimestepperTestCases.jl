#####
##### Loading and post-processing for the near-global case
#####
##### The near-global runs write three files per discretization -- `_dissipation`, `_surface` and `_average` --
##### plus the `_cost` record of `save_near_global_cost`. Contrarily to the idealized cases the fields are too
##### large to hold in memory for six variants at once, so `load_near_global` opens every time series with the
##### `OnDisk` backend and the diagnostics below reduce one snapshot at a time.
#####

using FFTW
using JLD2
using Oceananigans.Grids: φnodes, λnodes, znodes, xspacings
using Statistics: mean

near_global_path(folder, label, kind, prefix = "near_global") = joinpath(folder, "$(prefix)_$(label)_$(kind).jld2")

"""
    load_near_global(folder, label; architecture, prefix, dissipation, surface, average)

Open the output of the near-global run labelled `label` in `folder`.

$(SIGNATURES)

# Arguments
- `folder`: directory the run wrote into
- `label`: the discretization label of Table 1, one of `near_global_discretizations()`

# Keyword arguments
- `architecture`: architecture to load onto (default `CPU()`)
- `dissipation`, `surface`, `average`: whether to open each of the three output files, so that a notebook that
  only wants the surface maps does not pay for the buoyancy-variance budget

# Returns
- Dictionary with the buoyancy-variance dissipation `:Abx, :Aby, :Abz`, the squared gradients
  `:Gbx, :Gby, :Gbz`, the surface `:T, :S, :u, :v, :w`, the same fields time-averaged as `:T̄, :S̄, :ū, :v̄, :w̄`,
  the timings `:cost`, and the times of each group as `:times`, `:surface_times` and `:average_times`

Every field is opened with `backend = OnDisk()`: a single dissipation file is close to a gigabyte, so the
snapshots are read one at a time by the diagnostics rather than held. The dissipation and gradient fields are
already volume-integrated by the writer -- `Gbx = ∂x(b)² V` and likewise for the others -- which is why the
reductions below are plain sums and carry no metric of their own.
"""
function load_near_global(folder, label; architecture = CPU(), prefix = "near_global",
                          dissipation = true, surface = true, average = true)

    case = Dict{Symbol, Any}()
    case[:label] = label

    if dissipation
        path = near_global_path(folder, label, "dissipation", prefix)
        for name in (:Abx, :Aby, :Abz, :Gbx, :Gby, :Gbz)
            case[name] = FieldTimeSeries(path, string(name); architecture, backend = OnDisk())
        end
        case[:times] = case[:Abx].times
        case[:grid]  = case[:Abx].grid
    end

    if surface
        path = near_global_path(folder, label, "surface", prefix)
        for name in (:T, :S, :u, :v, :w)
            case[name] = FieldTimeSeries(path, string(name); architecture, backend = OnDisk())
        end
        case[:surface_times] = case[:u].times
    end

    if average
        path = near_global_path(folder, label, "average", prefix)
        for (name, averaged) in zip((:T, :S, :u, :v, :w), (:T̄, :S̄, :ū, :v̄, :w̄))
            case[averaged] = FieldTimeSeries(path, string(name); architecture, backend = OnDisk())
        end
        case[:average_times] = case[:ū].times
    end

    case[:cost] = load_near_global_cost(folder, label)

    return case
end

"""
    load_near_global_cases(folder; labels, prefix, kw...)

Open every near-global run present in `folder`, skipping the labels whose output is missing, so that a
notebook can be run while the remaining variants are still integrating.
"""
function load_near_global_cases(folder; labels = [d.label for d in near_global_discretizations()],
                                prefix = "near_global", kw...)
    available = filter(l -> isfile(near_global_path(folder, l, "dissipation", prefix)), labels)
    return available, [load_near_global(folder, l; prefix, kw...) for l in available]
end

#####
##### Timings
#####

"""
    load_near_global_cost(folder, label, prefix = "near_global")

The timings written by [`save_near_global_cost`](@ref), or `nothing` where the run predates the cost record.
"""
function load_near_global_cost(folder, label, prefix = "near_global")
    path = near_global_path(folder, label, "cost", prefix)
    isfile(path) || return nothing

    return JLD2.jldopen(path, "r") do file
        NamedTuple(Symbol(k) => file[k] for k in keys(file))
    end
end

"""
    near_global_cost_table(cases; reference = "AB2-SE")

Seconds per step, simulated years per wall-clock day (SYPD) at the production time step and the cost relative
to `reference`, for the cases whose timings are on disk. The relative column is the one Section 5.4.1 reads: a scheme that costs more per step
but takes a longer step can still finish the year first, which is what `simulated_years_per_day` shows.
"""
function near_global_cost_table(cases; reference = "AB2-SE")
    timed = filter(c -> !isnothing(c[:cost]), cases)
    isempty(timed) && return NamedTuple[]

    # a cost run steps below the time step its substeps are sized for, so the throughput is taken at the nominal one
    nominal_timestep(c) = get(c[:cost], :nominal_timestep, c[:cost].Δt)
    throughput(c) = near_global_simulated_years_per_day(nominal_timestep(c), c[:cost].seconds_per_step)

    reference_case = findfirst(c -> c[:label] == reference, timed)
    reference_throughput = isnothing(reference_case) ? NaN : throughput(timed[reference_case])

    return [(; label = c[:label],
               Δt = nominal_timestep(c),
               substeps = c[:cost].substeps,
               seconds_per_step = c[:cost].seconds_per_step,
               simulated_days_per_day = 365 * throughput(c),
               simulated_years_per_day = throughput(c),
               speedup = throughput(c) / reference_throughput) for c in timed]
end

#####
##### Numerical diffusivity
#####

"""
    near_global_diffusivity_profile(case, time_index = length(case[:times]))

The globally integrated numerical diffusivity `κ(z) = -Σ(Ax + Ay + Az) / 2Σ(Gx + Gy + Gz)`, with the sums taken
over the whole domain at each level, and the vertical component interpolated from faces to centers so that the
three directions are added at the same location.

$(SIGNATURES)

# Returns
- `(κ, z)`, the profile and the cell centers it sits on, both as plain vectors

Immersed cells contribute zero to both the numerator and the denominator, so the land mask needs no separate
treatment. Where a level is entirely land the ratio is `0/0` and comes back as `NaN`.
"""
function near_global_diffusivity_profile(case, time_index = length(case[:times]))
    Ax = horizontal_sum(case[:Abx][time_index])
    Ay = horizontal_sum(case[:Aby][time_index])
    Az = interpolate_to_center_z(horizontal_sum(case[:Abz][time_index]))

    Gx = horizontal_sum(case[:Gbx][time_index])
    Gy = horizontal_sum(case[:Gby][time_index])
    Gz = interpolate_to_center_z(horizontal_sum(case[:Gbz][time_index]))

    κ = @. - (Ax + Ay + Az) / (2 * (Gx + Gy + Gz))
    z = Array(znodes(case[:Abx].grid, Center()))

    return κ, z
end

"""
    near_global_diffusivity_map(case, time_index = length(case[:times]))

The depth-integrated numerical diffusivity as a map, `κ(λ, φ) = -Σ_z(Ax + Ay + Az) / 2Σ_z(Gx + Gy + Gz)`, with
the staggered directions interpolated onto cell centers. Returns `(κ, λ, φ)`.

The grid is curvilinear, so `λ` and `φ` come back as `Nx × Ny` matrices of the cell-center coordinates and not
as axes: a plot of `κ` reads them as the coordinates of each cell rather than as a rectangular mesh.
"""
function near_global_diffusivity_map(case, time_index = length(case[:times]))
    A = interpolate_to_center_x(vertical_sum(case[:Abx][time_index])) .+
        interpolate_to_center_y(vertical_sum(case[:Aby][time_index])) .+
        vertical_sum(case[:Abz][time_index])

    G = interpolate_to_center_x(vertical_sum(case[:Gbx][time_index])) .+
        interpolate_to_center_y(vertical_sum(case[:Gby][time_index])) .+
        vertical_sum(case[:Gbz][time_index])

    grid = case[:Abx].grid
    κ = @. - A / (2G)

    return κ, center_coordinates(grid)...
end

# λ and φ of the cell centers, as the `Nx × Ny` matrices a curvilinear grid carries
center_coordinates(grid) = Array(λnodes(grid, Center(), Center())), Array(φnodes(grid, Center(), Center()))

# Az and Δx on cell centers as plain matrices. The area weights and the zonal transform read these off the
# grid: on a curvilinear mesh neither of them factors into a function of latitude that could be rebuilt from
# the coordinates alone.
function center_areas(grid)
    Az = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Center),
                                                               Oceananigans.Operators.Az, grid)
    return Array(interior(compute!(Field(Az)), :, :, 1))
end

center_spacings(grid) = Array(interior(compute!(Field(xspacings(grid, Center(), Center(), Center()))), :, :, 1))

"""
    near_global_diffusivity_hovmoller(case; time_indices)

`κ(z, t)` for the whole output series, one column per output time, as the Hovmöller the other Results
notebooks draw. Costs one global reduction per output time, so pass a subsampled `time_indices` where the
series is long.
"""
function near_global_diffusivity_hovmoller(case; time_indices = eachindex(case[:times]))
    profiles = [near_global_diffusivity_profile(case, t)[1] for t in time_indices]
    z = Array(znodes(case[:Abx].grid, Center()))

    return reduce(hcat, profiles), z, case[:times][time_indices]
end

horizontal_sum(field) = Array(interior(compute!(Field(sum(field, dims = (1, 2)))), 1, 1, :))
vertical_sum(field)   = Array(interior(compute!(Field(sum(field, dims = 3))), :, :, 1))

interpolate_to_center_z(profile) = 0.5 .* (profile[1:end-1] .+ profile[2:end])
interpolate_to_center_y(a)       = 0.5 .* (a[:, 1:end-1] .+ a[:, 2:end])

# x is periodic on this configuration, so the face-to-center average wraps rather than losing a column
interpolate_to_center_x(a) = 0.5 .* (a .+ circshift(a, (-1, 0)))

#####
##### Surface diagnostics
#####

"""
    near_global_surface_speed(case, time_index = length(case[:surface_times]))

Surface speed `√(u² + v²)` on cell centers, with land as `NaN`. Returns `(speed, λ, φ)`, the coordinates as
the matrices of [`near_global_diffusivity_map`](@ref).

`u` and `v` are the grid-aligned components as the model carries them, so away from the Mercator part of the
mesh they are not eastward and northward; the speed is invariant under that rotation and needs no correction.
"""
function near_global_surface_speed(case, time_index = length(case[:surface_times]))
    u = interpolate_to_center_x(Array(interior(case[:u][time_index], :, :, 1)))
    v = interpolate_to_center_y_from_surface(case, time_index)

    grid = case[:u].grid
    speed = @. sqrt(u^2 + v^2)
    mask_land!(speed, case, time_index)

    return speed, center_coordinates(grid)...
end

"""
    near_global_surface_kinetic_energy(case, time_index = length(case[:surface_times]))

Surface kinetic energy `(u² + v²)/2` on cell centers, with land as `NaN`. Returns `(kinetic_energy, λ, φ)`.
"""
function near_global_surface_kinetic_energy(case, time_index = length(case[:surface_times]))
    speed, λ, φ = near_global_surface_speed(case, time_index)
    return speed .^ 2 ./ 2, λ, φ
end

"""
    near_global_eddy_kinetic_energy(case; time_indices)

Surface eddy kinetic energy, `⟨(u - ⟨u⟩)² + (v - ⟨v⟩)²⟩/2` with `⟨⟩` the average over `time_indices` of the
daily surface output. Returns `(eddy_kinetic_energy, λ, φ)`.

The average is accumulated snapshot by snapshot rather than by materializing the series, which for the daily
surface output of a two-year run is several gigabytes per variant.
"""
function near_global_eddy_kinetic_energy(case; time_indices = eachindex(case[:surface_times]))
    grid = case[:u].grid
    Nx, Ny, _ = size(grid)

    ū  = zeros(Nx, Ny); v̄  = zeros(Nx, Ny)
    u² = zeros(Nx, Ny); v² = zeros(Nx, Ny)

    for t in time_indices
        u = interpolate_to_center_x(Array(interior(case[:u][t], :, :, 1)))
        v = interpolate_to_center_y_from_surface(case, t)

        ū  .+= u;      v̄  .+= v
        u² .+= u .^ 2; v² .+= v .^ 2
    end

    n = length(time_indices)
    eddy_kinetic_energy = @. ((u² - ū^2 / n) + (v² - v̄^2 / n)) / (2n)
    mask_land!(eddy_kinetic_energy, case, first(time_indices))

    return eddy_kinetic_energy, center_coordinates(grid)...
end

"""
    near_global_surface_kinetic_energy_history(case; time_indices, latitude_band)

Area-weighted mean surface kinetic energy of every snapshot in `time_indices`, the time series the daily
surface output supports. Returns `(times, kinetic_energy)`.

$(SIGNATURES)

# Keyword arguments
- `time_indices`: the surface snapshots to read (default the whole series). One pass per case reads the daily
  output, so subsample where only the shape of the curve matters
- `latitude_band`: restrict the average to a band in degrees, `(-60, -56)` for the circumpolar band of
  [`near_global_zonal_spectrum`](@ref); `nothing` keeps the whole domain

The weights are the cell areas `Az` of the wet surface cells, read off the grid: on a curvilinear mesh the
area does not factor into a function of latitude, so there is nothing simpler than the metric itself to weight
with. The land mask is taken once, the bathymetry being fixed. This is the surface layer alone: the run writes
no three-dimensional velocity, so a volume-integrated kinetic energy is not available from the output on disk.
"""
function near_global_surface_kinetic_energy_history(case; time_indices = eachindex(case[:surface_times]),
                                                    latitude_band = nothing)

    grid = case[:u].grid
    φ = Array(φnodes(grid, Center(), Center()))
    Az = center_areas(grid)
    wet = Array(interior(case[:T][first(time_indices)], :, :, 1)) .!= 0

    in_band(i, j) = isnothing(latitude_band) || latitude_band[1] ≤ φ[i, j] ≤ latitude_band[2]
    weight = [wet[i, j] && in_band(i, j) ? Az[i, j] : 0.0 for i in axes(wet, 1), j in axes(wet, 2)]
    total_weight = sum(weight)

    kinetic_energy = map(time_indices) do t
        u = interpolate_to_center_x(Array(interior(case[:u][t], :, :, 1)))
        v = interpolate_to_center_y_from_surface(case, t)
        sum(@. weight * (u^2 + v^2) / 2) / total_weight
    end

    return case[:surface_times][time_indices], kinetic_energy
end

# v lives at (Center, Face, Center) and the meridional topology is face-extended, so the surface slice carries
# Ny+1 rows against the Ny of a centered field
function interpolate_to_center_y_from_surface(case, time_index)
    v = Array(interior(case[:v][time_index], :, :, 1))
    Ny = size(case[:u].grid, 2)
    return size(v, 2) == Ny + 1 ? interpolate_to_center_y(v) : v
end

# T is exactly zero only where the column is dry: the run floors the initial condition at -1.8°C and no wet
# surface cell reaches 0 to machine precision
function mask_land!(field, case, time_index)
    T = Array(interior(case[:T][time_index], :, :, 1))
    field[T .== 0] .= NaN
    return field
end

#####
##### Zonal spectra
#####

"""
    near_global_zonal_spectrum(case, name = :u; time_indices, latitude_band, minimum_wet_fraction)

Zonal power spectrum of a surface field, averaged over the rows of `latitude_band` and over the snapshots of
`time_indices`.

$(SIGNATURES)

# Arguments
- `case`: the dictionary of [`load_near_global`](@ref)
- `name`: which surface field, `:u`, `:v`, `:T` or `:S`

# Keyword arguments
- `time_indices`: the surface snapshots to average over, an index or a range (default the last snapshot). One
  snapshot carries a single realization per wavenumber, whose scatter is of the order of the estimate itself;
  averaging over `n` decorrelated snapshots brings it down as `1/√n`. The daily series is long, so the window
  is better taken over the whole stationary part of the run, subsampled past the eddy decorrelation time.
- `latitude_band`: the band to average over, in degrees. The default `(-60, -56)` is the band worth a zonal
  transform: it is periodic in the physical sense and entirely free of land -- south of Cape Horn and north of
  the Antarctic Peninsula -- so the transform sees an actual periodic signal rather than a coastline. Every
  latitude between 55°S and 45°S clips South America. It also sits in the interior of the domain, the nearest
  wall being the Antarctic coast some ten degrees further south, and in the Mercator part of the mesh, where
  a row is a line of constant latitude and the zonal direction is the grid direction.
- `minimum_wet_fraction`: rows with a smaller wet fraction are dropped. At the default of `1` only land-free
  rows are kept; relaxing it admits rows with a coastline, whose land cells are filled with the row's wet mean
  so that the transform sees a flat patch rather than a cliff to zero, and whose spectrum is contaminated at
  the short wavelengths accordingly.

# Returns
- `(spectrum, wavenumber)`, the power at each angular wavenumber in rad m⁻¹

Each row is transformed on its own spacing, read off the grid metric and contracting poleward, and the results
are interpolated onto the wavenumber axis of the band's central latitude before being averaged, so that the
poleward rows are not silently plotted against the equatorward axis.
"""
function near_global_zonal_spectrum(case, name = :u;
                                    time_indices = length(case[:surface_times]),
                                    latitude_band = (-60, -56),
                                    minimum_wet_fraction = 1)

    grid = case[:u].grid
    Nx = size(grid, 1)

    # One latitude and one spacing per row: within the band a row is a line of constant latitude along which
    # the spacing is constant, so the row means are those values and not an approximation of them
    φ  = vec(mean(Array(φnodes(grid, Center(), Center())), dims = 1))
    Δx = vec(mean(center_spacings(grid), dims = 1))

    land = Array(interior(case[:T][first(time_indices)], :, :, 1)) .== 0

    in_band(j) = latitude_band[1] ≤ φ[j] ≤ latitude_band[2]
    wet_fraction(j) = 1 - count(view(land, :, j)) / Nx

    rows = [j for j in eachindex(φ) if in_band(j) && wet_fraction(j) ≥ minimum_wet_fraction]

    if isempty(rows)
        wettest = maximum(wet_fraction, filter(in_band, eachindex(φ)); init = 0.0)
        throw(ArgumentError("No row in the band $latitude_band reaches a wet fraction of " *
                            "$minimum_wet_fraction; the wettest reaches $(round(wettest, digits = 3))"))
    end

    row_spacing(j) = Δx[j]

    spectrum, wavenumber = snapshot_zonal_spectrum(case, name, first(time_indices), rows, row_spacing)

    for time_index in Iterators.drop(time_indices, 1)
        snapshot, _ = snapshot_zonal_spectrum(case, name, time_index, rows, row_spacing)
        spectrum .+= snapshot
    end

    return spectrum ./ length(time_indices), wavenumber
end

function snapshot_zonal_spectrum(case, name, time_index, rows, row_spacing)
    data = name === :u ? interpolate_to_center_x(Array(interior(case[:u][time_index], :, :, 1))) :
           name === :v ? interpolate_to_center_y_from_surface(case, time_index) :
                         Array(interior(case[name][time_index], :, :, 1))

    land = Array(interior(case[:T][time_index], :, :, 1)) .== 0
    data = fill_land_with_row_mean(data, land, rows)

    reference_row = rows[cld(length(rows), 2)]
    reference_spectrum, reference_wavenumber = row_power_spectrum(data, reference_row, row_spacing(reference_row))

    spectrum = copy(reference_spectrum)
    for j in rows
        j == reference_row && continue
        row, wavenumber = row_power_spectrum(data, j, row_spacing(j))
        spectrum .+= interpolate_onto(row, wavenumber, reference_wavenumber)
    end

    return spectrum ./ length(rows), reference_wavenumber
end

function fill_land_with_row_mean(data, land, rows)
    filled = copy(data)

    for j in rows
        dry = view(land, :, j)
        any(dry) || continue
        filled[dry, j] .= sum(view(data, :, j)) / count(!, dry)
    end

    return filled
end

function row_power_spectrum(data, j, Δx)
    x = (0:size(data, 1)-1) .* Δx
    spectrum, wavenumber = power_spectrum_x(view(data, :, j), x)
    return real.(spectrum), wavenumber
end

function interpolate_onto(values, from, to)
    interpolated = similar(values, length(to))

    for (i, k) in enumerate(to)
        if k ≤ from[1]
            interpolated[i] = values[1]
        elseif k ≥ from[end]
            interpolated[i] = values[end]
        else
            m = searchsortedlast(from, k)
            w = (k - from[m]) / (from[m+1] - from[m])
            interpolated[i] = (1 - w) * values[m] + w * values[m+1]
        end
    end

    return interpolated
end
