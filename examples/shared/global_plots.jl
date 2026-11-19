function plot_global_cost(measurements, throughput)
    figure = Figure(size=(850, 430))
    axis   = Axis(figure[1, 1]; ylabel="Throughput relative to AB2-SE",
                xticks=(eachindex(measurements), [measurement.label for measurement in measurements]), xticklabelrotation=π/6)
    barplot!(axis, collect(eachindex(measurements)), throughput; color=[plot_style(measurement.label).color for measurement in measurements])
    hlines!(axis, [1]; color=:black, linestyle=:dash)
    return save_example(figure, "global_throughput")
end

const global_time_label = L"t\;[\mathrm{years}]"
const global_depth_label = L"z\;[\mathrm{km}]"
const global_diffusivity_label = L"\kappa\;[\mathrm{m^2\,s^{-1}}]"
const global_land_color = RGBf(0.82, 0.82, 0.82)

# Each cell is drawn between its four corners, unwrapped against the first one, so that a cell across the date
# line is not stretched over the map. The few cells wider than `maximum_cell_span` touch the pole and are dropped.
function tripolar_heatmap!(axis, values, λ, φ; longitude_limits=(0, 360), latitude_limits=(-90, 90), maximum_cell_span=90, kwargs...)
    vertices = Point2f[]
    colors   = Float32[]
    for j in axes(values, 2), i in axes(values, 1)
        isnan(values[i, j]) && continue
        first(latitude_limits) - 1 ≤ φ[i, j] ≤ last(latitude_limits) + 1 || continue
        corners   = ((λ[i, j], φ[i, j]), (λ[i+1, j], φ[i+1, j]), (λ[i+1, j+1], φ[i+1, j+1]), (λ[i, j+1], φ[i, j+1]))
        unwrapped = map(corner -> λ[i, j] + mod(first(corner) - λ[i, j] + 180, 360) - 180, corners)
        westmost, eastmost = extrema(unwrapped)
        eastmost - westmost > maximum_cell_span && continue
        for shift in (-360, 0, 360)
            eastmost + shift > first(longitude_limits) && westmost + shift < last(longitude_limits) || continue
            append!(vertices, Point2f(longitude + shift, last(corner)) for (longitude, corner) in zip(unwrapped, corners))
            append!(colors, ntuple(_ -> values[i, j], 4))
        end
    end
    faces = [face for q in 0:length(vertices)÷4-1 for face in (Makie.GLTriangleFace(4q + 1, 4q + 2, 4q + 3), Makie.GLTriangleFace(4q + 1, 4q + 3, 4q + 4))]
    return mesh!(axis, vertices, faces; color=colors, shading=NoShading, rasterize=2, kwargs...)
end

const global_wavenumber_label = L"k\;[\mathrm{rad\,m^{-1}}]"

function cubic_slope!(axis, k, energy)
    wavenumbers = k[end÷8:end÷2]
    slope       = 6 * energy[end÷8] * (wavenumbers / first(wavenumbers)) .^ (-3)
    lines!(axis, wavenumbers, slope; color=:gray40, linewidth=1)
    text!(axis, wavenumbers[end÷2], slope[end÷2]; text=L"k^{-3}", color=:gray40, align=(:left, :bottom), offset=(4, 4))
    return axis
end

positive_part(value) = value > 0 ? value : NaN
negative_part(value) = value < 0 ? -value : NaN

# κ < 0 is anti-diffusion: a logarithmic axis draws its magnitude as a separate dotted branch.
function signed_profile!(axis, κ, z; label=nothing, color, linestyle)
    lines!(axis, positive_part.(κ), z / 1000; color, linestyle, label)
    lines!(axis, negative_part.(κ), z / 1000; color, linestyle=:dot)
    return axis
end

function plot_global_diffusivity_profile(profiles)
    figure      = Figure(size=(1100, 560))
    logarithmic = (xscale=log10, xlabel=global_diffusivity_label)
    full_depth  = Axis(figure[1, 1]; title="Upper four kilometers", ylabel=global_depth_label, logarithmic...)
    upper_ocean = Axis(figure[1, 2]; title="Upper kilometer", xticks=([1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4],
                       [L"10^{-6}", L"2\cdot10^{-6}", L"5\cdot10^{-6}", L"10^{-5}", L"2\cdot10^{-5}", L"5\cdot10^{-5}", L"10^{-4}"]),
                       yticks=-1:0.2:0, logarithmic...)
    for profile in profiles, axis in (full_depth, upper_ocean)
        signed_profile!(axis, profile.κ, profile.z; label=publication_label(profile.label), plot_style(profile.label)...)
    end
    ylims!(full_depth, -4, 0)
    limits!(upper_ocean, 1e-6, 1e-4, -1, 0)
    Legend(figure[2, 1:2], full_depth; orientation=:horizontal)
    return save_example(figure, "global_diffusivity_profile")
end

function plot_global_diffusivity_wall_distance(profiles; diffusivity_limits=(1e-7, 3e-3))
    distances = axes(first(profiles).wall_diffusivity, 1)
    title(distance) = distance == last(distances) ? "More than $(distance - 1) cells from a wall" :
                      distance == 1 ? "1 cell from a wall" : "$distance cells from a wall"
    figure = Figure(size=(1100, 800))
    panels = map(distances) do distance
        row, column = fldmod1(distance, 3)
        axis = Axis(figure[row, column]; title=title(distance), xscale=log10, xticks=LogTicks(-6:-3), xlabel=row == 2 ? global_diffusivity_label : "",
                    xticklabelsvisible=row == 2, ylabel=column == 1 ? global_depth_label : "", yticklabelsvisible=column == 1)
        for profile in profiles
            signed_profile!(axis, view(profile.wall_diffusivity, distance, :), profile.z; label=publication_label(profile.label), plot_style(profile.label)...)
        end
        limits!(axis, diffusivity_limits..., -5, 0)
        return axis
    end
    Legend(figure[3, 1:3], first(panels); orientation=:horizontal)
    return save_example(figure, "global_diffusivity_wall_distance")
end

function plot_global_diffusivity_section(sections; diffusivity_limits=(1e-6, 1e-3), difference_limit=1e-3)
    reference         = only(section for section in sections if section.label == "WRK3-SE")
    differences       = filter(section -> section.label != reference.label, sections)
    figure            = Figure(size=(1100, 760))
    positive_levels   = 10.0 .^ range(log10(difference_limit) - 4, log10(difference_limit); length=13)
    difference_levels = vcat(-reverse(positive_levels), positive_levels)
    difference_scale  = Makie.Symlog10(first(positive_levels); linscale=10/9)

    for (index, section) in enumerate((reference, differences...))
        row, column = fldmod1(index, 3)
        is_reference = index == 1
        title        = publication_label(section.label) * (is_reference ? "" : " − WRK3-SE")
        axis         = Axis(figure[row, column]; title, backgroundcolor=global_land_color, xgridvisible=false, ygridvisible=false,
                    xlabel=row == 2 ? L"\varphi\;[^\circ\mathrm{N}]" : "", xticks=-60:30:60, xticklabelsvisible=row == 2,
                    ylabel=column == 1 ? global_depth_label : "", yticklabelsvisible=column == 1)
        if is_reference
            levels = 10.0 .^ range(log10(first(diffusivity_limits)), log10(last(diffusivity_limits)); length=17)
            values = max.(section.κ, first(diffusivity_limits) / 10)
            scale, colormap = log10, :YlOrRd
        else
            @assert section.latitudes == reference.latitudes && section.z == reference.z
            levels = difference_levels
            values = section.κ - reference.κ
            scale, colormap = difference_scale, :bwr
        end
        # Transform values and levels together so the colorbar preserves the white band around zero.
        bands = contourf!(axis, section.latitudes, section.z / 1000, scale.(values); levels=scale.(levels), colormap,
                          extendlow=:auto, extendhigh=:auto)
        contour!(axis, section.latitudes, section.z / 1000, scale.(values); levels=scale.(levels), color=(:white, 0.4), linewidth=0.25)
        limits!(axis, -80, 90, -5, 0)
        if is_reference
            ticks = (log10.([1e-6, 1e-5, 1e-4, 1e-3]), [L"10^{-6}", L"10^{-5}", L"10^{-4}", L"10^{-3}"])
            Colorbar(figure[3, 1], bands; vertical=false, ticks, label=global_diffusivity_label, height=12)
        elseif index == length(sections)
            powers      = 10.0 .^ (-7:2:-3)
            ticks       = vcat(-reverse(powers), 0, powers)
            tick_labels = signed_power_label.(ticks)
            Colorbar(figure[3, 2:3], bands; vertical=false, ticks=(scale.(ticks), tick_labels), height=12,
                     label=L"\kappa-\kappa^{\mathrm{WRK3-SE}}\;[\mathrm{m^2\,s^{-1}}]")
        end
    end
    colgap!(figure.layout, 18)
    rowgap!(figure.layout, 18)
    return save_example(figure, "global_diffusivity_section")
end

# The windows are shown at their midpoints, and a line stops where a window is missing from the record.
function plot_global_diffusivity_history(histories)
    layers = ((:κ, "Whole domain"), (:deep_diffusivity, "Below one kilometer"))
    figure = Figure(size=(1100, 400))
    panels = map(enumerate(layers)) do (column, (name, title))
        axis = Axis(figure[1, column]; title, xlabel=global_time_label, ylabel=column == 1 ? L"\kappa\;[10^{-5}\,\mathrm{m^2\,s^{-1}}]" : "")
        for history in histories
            style = plot_style(history.label)
            scatterlines!(axis, history.years, 1e5 * getproperty(history, name); style..., markercolor=style.color, markersize=9, label=publication_label(history.label))
        end
        return axis
    end
    Legend(figure[2, 1:2], first(panels); orientation=:horizontal)
    return save_example(figure, "global_diffusivity_history")
end

function plot_global_surface_tracer_histories(histories)
    figure                 = Figure(size=(1000, 640))
    temperature            = Axis(figure[1, 1]; title="Surface temperature", ylabel=L"\langle T\rangle_A\;[^\circ\mathrm{C}]", xticklabelsvisible=false)
    salinity               = Axis(figure[1, 2]; title="Surface salinity", ylabel=L"\langle S\rangle_A\;[\mathrm{psu}]", xticklabelsvisible=false)
    temperature_difference = Axis(figure[2, 1]; title="Difference from AB2-SE", xlabel=global_time_label,
                                  ylabel=L"\Delta\langle T\rangle_A\;[10^{-3}\,{}^\circ\mathrm{C}]")
    salinity_difference = Axis(figure[2, 2]; title="Difference from AB2-SE", xlabel=global_time_label,
                               ylabel=L"\Delta\langle S\rangle_A\;[10^{-3}\,\mathrm{psu}]")
    for history in histories
        style = plot_style(history.label)
        label = publication_label(history.label)
        lines!(temperature, history.years, history.temperature; style..., label)
        lines!(salinity, history.years, history.salinity; style...)
        lines!(temperature_difference, history.years, 1e3 * history.temperature_difference; style...)
        lines!(salinity_difference, history.years, 1e3 * history.salinity_difference; style...)
    end
    Legend(figure[3, 1:2], temperature; orientation=:horizontal)
    return save_example(figure, "global_surface_tracer_histories")
end

function plot_global_surface(surfaces, regions)
    figure = Figure(size=(1200, 1120), fontsize=22)
    Colorbar(figure[1, 1]; colormap=:magma, limits=(0, 0.5), vertical=false, height=9,
             width=Relative(0.8), flipaxis=true, ticks=0:0.1:0.5, ticklabelsize=18,
             label=L"\sqrt{u^2+v^2}\;[\mathrm{m\,s^{-1}}]")
    for (row, label) in enumerate(("AB2-SE", "WRK3-SE", "MRK4-SE"))
        surface = only(surface for surface in surfaces if surface.label == label)
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
                   label=publication_label(surface.label))
        end
        spectrum = first(surfaces).spectra[row]
        cubic_slope!(axis, spectrum.wavenumbers, vec(mean(spectrum.energy; dims=2)))
        row == 1 && Legend(figure[5, 1:2], axis; orientation=:horizontal, nbanks=2, framevisible=true)
    end
    colsize!(figure.layout, 1, Relative(0.65))
    colgap!(figure.layout, 28)
    rowgap!(figure.layout, 18)
    return save_example(figure, "global_surface_comparison")
end
