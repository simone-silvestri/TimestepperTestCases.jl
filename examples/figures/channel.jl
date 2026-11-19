# # Equilibrated channel: numerical diffusivity and dissipation
#
# The input files contain the averaged buoyancy-variance tendencies and squared gradients. This example uses
# the final saved averaging window, recorded below for each case. Only the averages files are required;
# three-dimensional snapshots and restart files are unnecessary for these diagnostics.

using TimestepperTestCases
using Oceananigans
using Statistics
include("../shared/data_sources.jl") #src
include("../shared/compact_data.jl") #src
include("../shared/plotting.jl") #src
function plot_channel_diffusivity(profiles)
    figure = Figure(size=(595, 210), fontsize=10)
    labels = (L"|\kappa_{\mathrm{num},h}|\;[\mathrm{m^2\,s^{-1}}]", L"|\kappa_{\mathrm{num}}|\;[\mathrm{m^2\,s^{-1}}]",
              L"\kappa_{\mathrm{phys},z}/\kappa_{\mathrm{num}}")
    axes = [Axis(figure[1, column]; title=label, ylabel=column == 1 ? L"z\;[\mathrm{km}]" : "", xscale=log10,
                 yticklabelsvisible=column == 1, yticks=-3:0.5:0, limits=(nothing, nothing, -3, 0),
                 xgridcolor=(:black, 0.08), ygridcolor=(:black, 0.08)) for (column, label) in enumerate(labels)]
    axes[1].xticks = ([1, 10], [L"10^0", L"10^1"])
    axes[2].xticks = ([1e-6, 1e-5, 1e-4, 1e-3], [L"10^{-6}", L"10^{-5}", L"10^{-4}", L"10^{-3}"])
    axes[3].xticks = ([0.1, 1, 10], [L"10^{-1}", L"10^0", L"10^1"])
    xlims!(axes[2], 5e-7, 1e-3)
    for profile in profiles
        style = plot_style(profile.label; distinct_colors=true)
        label = publication_label(profile.label)
        for (column, values) in enumerate((profile.κₕ, profile.κ, profile.ratio))
            displayed = column < 3 ? abs.(values) : [value > 0 ? value : NaN for value in values]
            lines!(axes[column], displayed, profile.z / 1000; style..., linewidth=1.4, label)
        end
    end
    vlines!(axes[3], [1]; color=(:black, 0.35), linestyle=:dash, linewidth=0.8)
    Legend(figure[1, 4], axes[1]; framevisible=true, patchsize=(20, 8), rowgap=3, padding=(5, 5, 5, 5))
    colgap!(figure.layout, 12)
    return save_example(figure, "channel_diffusivity"; pt_per_unit=1)
end

function plot_channel_dissipation(reference, differences)
    figure = Figure(size=(764, 383), fontsize=10)
    for (row, property, difference_property) in ((1, :Pₕ, :ΔPₕ), (2, :Pz, :ΔPz))
        for (column, section) in enumerate((reference, differences...))
            is_reference = column == 1
            values       = getproperty(section, is_reference ? property : difference_property)
            lower, upper, count = is_reference ? (-9.5, -4.5, 14) : row == 1 ? (-7.5, -4.5, 13) : (-8.5, -5.5, 13)
            positive_levels = 10.0 .^ range(lower, upper; length=count)
            levels          = vcat(-reverse(positive_levels), 0, positive_levels)
            scale           = Makie.Symlog10(first(positive_levels); linscale=10/9)
            title           = row == 1 ? (is_reference ? section.label : section.label * " − WRK3-SE") : ""
            axis            = Axis(figure[row, is_reference ? 1 : column+1]; title, xlabel=row == 2 ? L"y\;(\mathrm{km})" : "",
                        ylabel=is_reference ? "Depth (km)" : "", yticklabelsvisible=is_reference, xticklabelsvisible=row == 2,
                        yticks=-3:1:0, xticks=500:500:1500, limits=(0, 2000, -3, 0), xgridvisible=false, ygridvisible=false)
            bands = contourf!(axis, section.y / 1000, section.z / 1000, values; levels, colorscale=scale,
                              colormap=:bwr, extendlow=:auto, extendhigh=:auto)
            contour!(axis, section.y / 1000, section.z / 1000, values; levels, color=(:white, 0.4), linewidth=0.25)
            if is_reference || column == 4
                powers      = 10.0 .^ (ceil(Int, lower):floor(Int, upper))
                ticks       = vcat(-reverse(powers), 0, powers)
                tick_labels = signed_power_label.(ticks)
                label       = is_reference ? (row == 1 ? L"\overline{P_x+P_y}" : L"\overline{P_z}") :
                                       (row == 1 ? L"\Delta\overline{P_x+P_y}" : L"\Delta\overline{P_z}")
                Colorbar(figure[row, is_reference ? 2 : 6], bands; scale, ticks=(ticks, tick_labels), label, width=6,
                         ticklabelsize=7, ticksize=2, labelsize=9)
            end
        end
    end
    colgap!(figure.layout, 8)
    rowgap!(figure.layout, 8)
    return save_example(figure, "channel_dissipation"; pt_per_unit=1)
end

directory = isempty(ARGS) ? nothing : first(ARGS)
labels = ("AB2-SE", "WRK3-SE-SM05", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM", "WRK3-UP")
profiles = NamedTuple[]
sections = NamedTuple[]

# The flux components live on different faces. Interpolation to cell centers precedes their sum. Numerical
# diffusivity is the ratio of the horizontally averaged variance tendency to twice the averaged squared gradient.

for label in labels
    data = compact_dataset("channel_simulation", "averages_$label.jld2", prepare_channel; directory)
    @info "Channel averaging window" case=label final_year=data.time / (360 * 86400)
    Pₕ = data.fields["Abx"] + data.fields["Aby"]
    Pz = data.fields["Abz"]
    Gₕ = data.fields["Gbx"] + data.fields["Gby"]
    Gz = data.fields["Gbz"]
    κₕ = vec(-sum(Pₕ; dims=1) ./ (2sum(Gₕ; dims=1)))
    κ  = vec(-sum(Pₕ + Pz; dims=1) ./ (2sum(Gₕ + Gz; dims=1)))
    κᵖ = vec(-sum(data.fields["Dbz"]; dims=1) ./ (2sum(Gz; dims=1)))
    push!(profiles, (; label, data.z, κₕ, κ, ratio=κᵖ ./ κ, data.time))
    push!(sections, (; label, data.y, data.z, Pₕ, Pz))
end

# The section differences retain the sign of the variance tendency and use WRK3-SE as the reference.
reference = only(filter(section -> section.label == "WRK3-SE", sections))
differences = [(; section..., ΔPₕ=section.Pₕ - reference.Pₕ, ΔPz=section.Pz - reference.Pz) for section in sections if section.label in ("AB2-SE", "MRK4-SE", "WRK3-UP")]

plot_channel_diffusivity(profiles) #src
plot_channel_dissipation(reference, differences) #src
# ![Channel numerical diffusivity](channel_diffusivity.png)
# ![Channel variance-tendency differences](channel_dissipation.png)
