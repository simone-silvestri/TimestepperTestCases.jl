# # Global ocean: surface speed and regional velocity spectra
#
# Run `julia --project=examples examples/global_surface_results.jl global_ocean/global` from the repository root.
# The maps use the final common daily snapshot. Velocity spectra use 96-cell, land-free sections, with a linear
# detrend and Hann window, sampled every five days over the final 365 days. Each region averages the two velocity spectra.

using TimestepperTestCases
using Oceananigans
using JLD2
using Statistics
include("compact_data.jl")
include("plotting.jl")

directory = isempty(ARGS) ? joinpath(@__DIR__, "..", "global_ocean", "global") : first(ARGS)
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
figure = Figure(size=(1200, 1120), fontsize=22)
Colorbar(figure[1, 1]; colormap=:magma, limits=(0, 0.5), vertical=false, height=9,
         width=Relative(0.8), flipaxis=true, ticks=0:0.1:0.5, ticklabelsize=18,
         label=L"\sqrt{u^2+v^2}\;[\mathrm{m\,s^{-1}}]")
for (row, index) in enumerate((1, 3, 5))
    surface = surfaces[index]
    axis    = Axis(figure[row+1, 1]; title="($(('a', 'b', 'c')[row])) $(surface.label)", titlealign=:left,
                backgroundcolor=global_land_color, xgridvisible=false, ygridvisible=false,
                xlabel=row == 3 ? L"\lambda\;[^\circ\mathrm{E}]" : "", ylabel=L"\varphi\;[^\circ\mathrm{N}]",
                xticks=0:90:360, yticks=-60:30:60, xticklabelsvisible=row == 3)
    tripolar_heatmap!(axis, sqrt.(2surface.kinetic_energy), surface.λ, surface.φ; colormap=:magma, colorrange=(0, 0.5))
    limits!(axis, 0, 360, -80, 85)
    for (region, abbreviation) in zip(regions, ("GS", "KE", "SO"))
        west, east = region.longitudes
        south, north = region.latitudes
        lines!(axis, [west, east, east, west, west], [south, south, north, north, south]; color=:cyan, linewidth=1.7)
        text!(axis, (west+east)/2, north; text=abbreviation, color=:cyan, align=(:center, :bottom), fontsize=17)
    end
end
for (row, region) in enumerate(regions)
    axis = Axis(figure[row+1, 2]; title="($(('d', 'e', 'f')[row])) $(region.name)", titlealign=:left,
                xlabel=row == 3 ? global_wavenumber_label : "", xscale=log10, yscale=log10,
                ylabel=L"E_k(k)\;[\mathrm{m^3\,s^{-2}}]", xticks=LogTicks(-5:-4),
                xticklabelsvisible=row == 3, xgridcolor=(:black, 0.07), ygridcolor=(:black, 0.07))
    for surface in surfaces
        spectrum = surface.spectra[row]
        lines!(axis, spectrum.wavenumbers[3:end], vec(mean(spectrum.energy; dims=2))[3:end];
               plot_style(surface.label; distinct_colors=true)..., linewidth=2,
               label=replace(surface.label, "SM05" => "F1"))
    end
    spectrum = first(surfaces).spectra[row]
    cubic_slope!(axis, spectrum.wavenumbers, vec(mean(spectrum.energy; dims=2)))
    row == 1 && Legend(figure[5, 1:2], axis; orientation=:horizontal, nbanks=2, framevisible=true)
end
colsize!(figure.layout, 1, Relative(0.65))
colgap!(figure.layout, 28)
rowgap!(figure.layout, 18)
save_example(figure, "global_surface_comparison")
