# # Global ocean: surface speed and regional velocity spectra
#
# Run `julia --project=examples examples/figures/global_surface_results.jl /path/to/global_ocean` from the repository root.
# The maps use the final common daily snapshot. Velocity spectra use 96-cell, land-free sections, with a linear
# detrend and Hann window, sampled every five days over the final 365 days. Each region averages the two velocity spectra.

using TimestepperTestCases
using Oceananigans
using JLD2
using Statistics
include("../shared/compact_data.jl") #src
include("../shared/plotting.jl") #src
include("../shared/global_plots.jl") #src

directory = isempty(ARGS) ? joinpath(@__DIR__, "..", "..", "global_ocean") : first(ARGS)
labels = ["AB2-SE", "WRK3-SE-SM05", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-UP"]
regions = [(; name="Gulf Stream", longitudes=(280, 320), latitudes=(28, 50)),
           (; name="Kuroshio Extension", longitudes=(135, 175), latitudes=(25, 47)),
           (; name="Southern Ocean", longitudes=(30, 90), latitudes=(-70, -38))]

surfaces = map(labels) do label
    path     = joinpath(directory, "global_ocean_$(label)_surface.jld2")
    diagnostics = joinpath(directory, "global_ocean_$(label)_velocity_diagnostics.jld2")
    if isfile(diagnostics) && mtime(diagnostics) ≥ mtime(path)
        @info "Reading surface velocity diagnostics" label
        return JLD2.load(diagnostics, "result")
    end
    geometry = global_ocean_geometry(path)
    @info "Reducing surface velocities" label
    result = jldopen(path) do file
        iterations, days = daily_records(file)
        selected       = findall(≥(last(days) - 365), days)[1:5:end]
        last_iteration = last(iterations)
        kinetic_energy = surface_kinetic_energy(file["timeseries/u/$last_iteration"], file["timeseries/v/$last_iteration"], geometry)
        spectra        = map(regions) do region
            segments = land_free_segments(region, geometry; segment_length=96)
            spacings = [mean(geometry.Δx[i, segment.row] for i in segment.columns) for segment in segments]
            spectrum(field) = mean(tapered_spectrum(Float64[field[geometry.Hx+i, geometry.Hy+segment.row, 1] for i in segment.columns], spacing)
                                   for (segment, spacing) in zip(segments, spacings))
            energy = stack((spectrum(file["timeseries/u/$(iterations[n])"]) + spectrum(file["timeseries/v/$(iterations[n])"])) / 2 for n in selected)
            return (; region.name, wavenumbers=2π * (0:47) / (96mean(spacings)), energy, number_of_segments=length(segments))
        end
        return (; label, snapshot_day=last(days), spectral_days=days[selected], kinetic_energy, spectra, λ=geometry.λᶠᶠ, φ=geometry.φᶠᶠ)
    end
    @info "Surface reduction complete" label result.snapshot_day length(result.spectral_days)
    jldsave(diagnostics; result)
    return result
end

@assert all(surface -> surface.snapshot_day == first(surfaces).snapshot_day && surface.spectral_days ≈ first(surfaces).spectral_days, surfaces)
plot_global_surface(surfaces, regions) #src
# ![Surface speed and regional velocity spectra](global_surface_comparison.png)
