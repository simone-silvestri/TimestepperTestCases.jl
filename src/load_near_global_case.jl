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
using Oceananigans.Grids: φnodes, λnodes, znodes

near_global_dissipation_path(folder, label) = joinpath(folder, "near_global_$(label)_dissipation.jld2")
near_global_surface_path(folder, label)     = joinpath(folder, "near_global_$(label)_surface.jld2")
near_global_average_path(folder, label)     = joinpath(folder, "near_global_$(label)_average.jld2")
near_global_cost_path(folder, label)        = joinpath(folder, "near_global_$(label)_cost.jld2")

"""
    load_near_global(folder, label; architecture, dissipation, surface, average)

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
function load_near_global(folder, label; architecture = CPU(),
                          dissipation = true, surface = true, average = true)

    case = Dict{Symbol, Any}()
    case[:label] = label

    if dissipation
        path = near_global_dissipation_path(folder, label)
        for name in (:Abx, :Aby, :Abz, :Gbx, :Gby, :Gbz)
            case[name] = FieldTimeSeries(path, string(name); architecture, backend = OnDisk())
        end
        case[:times] = case[:Abx].times
        case[:grid]  = case[:Abx].grid
    end

    if surface
        path = near_global_surface_path(folder, label)
        for name in (:T, :S, :u, :v, :w)
            case[name] = FieldTimeSeries(path, string(name); architecture, backend = OnDisk())
        end
        case[:surface_times] = case[:u].times
    end

    if average
        path = near_global_average_path(folder, label)
        for (name, averaged) in zip((:T, :S, :u, :v, :w), (:T̄, :S̄, :ū, :v̄, :w̄))
            case[averaged] = FieldTimeSeries(path, string(name); architecture, backend = OnDisk())
        end
        case[:average_times] = case[:ū].times
    end

    case[:cost] = load_near_global_cost(folder, label)

    return case
end

"""
    load_near_global_cases(folder; labels, kw...)

Open every near-global run present in `folder`, skipping the labels whose output is missing, so that a
notebook can be run while the remaining variants are still integrating.
"""
function load_near_global_cases(folder; labels = [d.label for d in near_global_discretizations()], kw...)
    available = filter(l -> isfile(near_global_dissipation_path(folder, l)), labels)
    return available, [load_near_global(folder, l; kw...) for l in available]
end

#####
##### Timings
#####

"""
    load_near_global_cost(folder, label)

The timings written by [`save_near_global_cost`](@ref), or `nothing` where the run predates the cost record.
"""
function load_near_global_cost(folder, label)
    path = near_global_cost_path(folder, label)
    isfile(path) || return nothing

    return JLD2.jldopen(path, "r") do file
        NamedTuple(Symbol(k) => file[k] for k in keys(file))
    end
end

"""
    near_global_cost_table(cases; reference = "AB2-SE")

Seconds per step, simulated days per wall-clock day and the cost relative to `reference`, for the cases whose
timings are on disk. The relative column is the one Section 5.4.1 reads: a scheme that costs more per step
but takes a longer step can still finish the year first, which is what `simulated_days_per_day` shows.
"""
function near_global_cost_table(cases; reference = "AB2-SE")
    timed = filter(c -> !isnothing(c[:cost]), cases)
    isempty(timed) && return NamedTuple[]

    # simulated seconds per wall second, which is the same number as simulated days per wall-clock day
    throughput(c) = c[:cost].Δt / c[:cost].seconds_per_step

    reference_case = findfirst(c -> c[:label] == reference, timed)
    reference_throughput = isnothing(reference_case) ? NaN : throughput(timed[reference_case])

    return [(; label = c[:label],
               Δt = c[:cost].Δt,
               substeps = c[:cost].substeps,
               seconds_per_step = c[:cost].seconds_per_step,
               simulated_days_per_day = throughput(c),
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

    return κ, Array(λnodes(grid, Center())), Array(φnodes(grid, Center()))
end

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

Surface speed `√(u² + v²)` on cell centers, with land as `NaN`. Returns `(speed, λ, φ)`.
"""
function near_global_surface_speed(case, time_index = length(case[:surface_times]))
    u = interpolate_to_center_x(Array(interior(case[:u][time_index], :, :, 1)))
    v = interpolate_to_center_y_from_surface(case, time_index)

    grid = case[:u].grid
    speed = @. sqrt(u^2 + v^2)
    mask_land!(speed, case, time_index)

    return speed, Array(λnodes(grid, Center())), Array(φnodes(grid, Center()))
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

    return eddy_kinetic_energy, Array(λnodes(grid, Center())), Array(φnodes(grid, Center()))
end

# v lives at (Center, Face, Center) and the meridional direction is bounded, so the surface slice carries Ny+1
# rows against the Ny of a centered field
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
    near_global_zonal_spectrum(case, name = :u; time_index, latitude_band, skip_land)

Zonal power spectrum of a surface field, averaged over the rows of `latitude_band`.

$(SIGNATURES)

# Arguments
- `case`: the dictionary of [`load_near_global`](@ref)
- `name`: which surface field, `:u`, `:v`, `:T` or `:S`

# Keyword arguments
- `time_index`: which surface snapshot (default the last)
- `latitude_band`: the band to average over, in degrees. The default `(-60, -56)` is the one band of this
  configuration that is worth a zonal transform: it is periodic in the physical sense and, uniquely, entirely
  free of land -- south of Cape Horn and north of the Antarctic Peninsula -- so the transform sees an actual
  periodic signal rather than a coastline. Every latitude between 55°S and 45°S clips South America.
- `minimum_wet_fraction`: rows with a smaller wet fraction are dropped. At the default of `1` only land-free
  rows are kept; relaxing it admits rows with a coastline, whose land cells are filled with the row's wet mean
  so that the transform sees a flat patch rather than a cliff to zero, and whose spectrum is contaminated at
  the short wavelengths accordingly.

# Returns
- `(spectrum, wavenumber)`, the power at each angular wavenumber in rad m⁻¹

Each row is transformed on its own spacing -- `Δx = 2πR cos φ / Nx` shrinks poleward -- and the results are
interpolated onto the wavenumber axis of the band's central latitude before being averaged, so that the
poleward rows are not silently plotted against the equatorward axis.
"""
function near_global_zonal_spectrum(case, name = :u;
                                    time_index = length(case[:surface_times]),
                                    latitude_band = (-60, -56),
                                    minimum_wet_fraction = 1)

    grid = case[:u].grid
    φ = Array(φnodes(grid, Center()))
    Nx = size(grid, 1)

    data = name === :u ? interpolate_to_center_x(Array(interior(case[:u][time_index], :, :, 1))) :
           name === :v ? interpolate_to_center_y_from_surface(case, time_index) :
                         Array(interior(case[name][time_index], :, :, 1))

    land = Array(interior(case[:T][time_index], :, :, 1)) .== 0

    in_band(j) = latitude_band[1] ≤ φ[j] ≤ latitude_band[2]
    wet_fraction(j) = 1 - count(view(land, :, j)) / Nx

    rows = [j for j in eachindex(φ) if in_band(j) && wet_fraction(j) ≥ minimum_wet_fraction]

    if isempty(rows)
        wettest = maximum(wet_fraction, filter(in_band, eachindex(φ)); init = 0.0)
        throw(ArgumentError("No row in the band $latitude_band reaches a wet fraction of " *
                            "$minimum_wet_fraction; the wettest reaches $(round(wettest, digits = 3))"))
    end

    data = fill_land_with_row_mean(data, land, rows)

    radius = grid_radius(grid)
    row_spacing(j) = 2π * radius * cosd(φ[j]) / Nx

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

grid_radius(grid) = hasproperty(grid, :radius) ? grid.radius : grid_radius(grid.underlying_grid)

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
