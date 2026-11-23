# # Global ocean: surface speed and regional velocity spectra
#
# In this example, we draw the surface figure of the manuscript (Figure 17) from the six `global_ocean_<label>_surface.jld2`
# files of the record, which hold the daily surface velocities of the final 365 days. Run it from the repository root
# with `julia --project=examples examples/figures/global_surface_results.jl /path/to/global_ocean`.

using TimestepperTestCases
using Oceananigans
using JLD2
using Statistics
include("../shared/compact_data.jl") #src
include("../shared/plotting.jl") #src
include("../shared/global_plots.jl") #src

directory = isempty(ARGS) ? joinpath(@__DIR__, "..", "..", "global_ocean") : first(ARGS)
labels = ["AB2-SE", "WRK3-SE-SM05", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-UP"]

# ## Regions
#
# Numerical dissipation removes energy first at the smallest resolved scales, so kinetic-energy spectra show where each
# scheme keeps its energy. We compute them in three eddy-rich regions, drawn as rectangles on the maps.

regions = [(; name="Gulf Stream", longitudes=(280, 320), latitudes=(28, 50)),
           (; name="Kuroshio Extension", longitudes=(135, 175), latitudes=(25, 47)),
           (; name="Southern Ocean", longitudes=(30, 90), latitudes=(-70, -38))]

# ## Spectrum of a region
#
# Inside a region, the spectra are taken along grid rows, over land-free segments of 96 cells that overlap by three
# quarters. A cell enters a segment only if the velocity faces on its western and southern sides separate two wet cells.
# Each segment is linearly detrended and multiplied by a Hann window before its spectral density is computed. The
# spectrum of a snapshot averages the segments and the two velocity components; `energy` holds one column per snapshot.

function region_spectrum(file, iterations, selected, region, geometry)
    segments = land_free_segments(region, geometry; segment_length=96)
    spacings = [mean(geometry.Δx[i, segment.row] for i in segment.columns) for segment in segments]
    spectrum(field) = mean(tapered_spectrum(Float64[field[geometry.Hx+i, geometry.Hy+segment.row, 1] for i in segment.columns], spacing)
                           for (segment, spacing) in zip(segments, spacings))
    energy = stack((spectrum(file["timeseries/u/$(iterations[n])"]) + spectrum(file["timeseries/v/$(iterations[n])"])) / 2 for n in selected)
    return (; region.name, wavenumbers=2π * (0:47) / (96mean(spacings)), energy, number_of_segments=length(segments))
end

# ## Reducing one configuration
#
# Restarts leave extra records shortly after a whole day, so `daily_records` keeps the record closest to each day. The
# maps use the surface kinetic energy of the final snapshot, and the spectra every fifth daily snapshot of the final
# 365 days. The reduction is cached next to the input as `global_ocean_<label>_velocity_diagnostics.jld2` and reused
# while it is newer than the surface file.

function reduce_surface(label)
    path        = joinpath(directory, "global_ocean_$(label)_surface.jld2")
    diagnostics = joinpath(directory, "global_ocean_$(label)_velocity_diagnostics.jld2")
    if isfile(diagnostics) && mtime(diagnostics) ≥ mtime(path)
        @info "Reading surface velocity diagnostics" label
        return JLD2.load(diagnostics, "result")
    end
    geometry = global_ocean_geometry(path)
    @info "Reducing surface velocities" label
    result = jldopen(path) do file
        iterations, days = daily_records(file)
        selected         = findall(≥(last(days) - 365), days)[1:5:end]
        last_iteration   = last(iterations)
        kinetic_energy   = surface_kinetic_energy(file["timeseries/u/$last_iteration"], file["timeseries/v/$last_iteration"], geometry)
        spectra          = map(region -> region_spectrum(file, iterations, selected, region, geometry), regions)
        return (; label, snapshot_day=last(days), spectral_days=days[selected], kinetic_energy, spectra, λ=geometry.λᶠᶠ, φ=geometry.φᶠᶠ)
    end
    @info "Surface reduction complete" label result.snapshot_day length(result.spectral_days)
    jldsave(diagnostics; result)
    return result
end

surfaces = map(reduce_surface, labels)

# Every configuration is compared on the same snapshot day and the same spectral days.

@assert all(surface -> surface.snapshot_day == first(surfaces).snapshot_day && surface.spectral_days ≈ first(surfaces).spectral_days, surfaces)

plot_global_surface(surfaces, regions) #src
# ![Surface speed and regional velocity spectra](global_surface_comparison.png)
#
# All cases show the major currents and a similar eddy field, and their spectra are close in all three regions. WRK3-UP
# has less energy at every wavenumber, consistent with its larger numerical diffusivity, while the Runge–Kutta schemes
# show no systematic gain of surface kinetic energy over AB2-SE.
