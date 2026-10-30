using CairoMakie
using LaTeXStrings

set_theme!(Theme(fontsize=16, fonts=(; regular="TeX Gyre Termes"), Legend=(; framevisible=true)))

function save_example(figure, name)
    directory = joinpath(@__DIR__, "output")
    mkpath(directory)
    for extension in ("png", "pdf")
        save(joinpath(directory, name * "." * extension), figure)
    end
    return figure
end

function plot_internal_tide(profiles, differences)
    figure = Figure(size=(1100, 720))
    history = Axis(figure[1, 1:2]; xlabel=L"t\;[\mathrm{days}]", ylabel=L"\langle\kappa\rangle\;[\mathrm{m^2\,s^{-1}}]")
    profile_axis = Axis(figure[1, 3]; xlabel=L"\overline{\kappa}\;[\mathrm{m^2\,s^{-1}}]", ylabel=L"z\;[\mathrm{km}]")
    for profile in profiles
        style = plot_style(profile.label)
        lines!(history, profile.t / 86400, profile.κᶻ; style..., label=profile.label)
        lines!(profile_axis, profile.κᵗ, profile.z / 1000; style...)
    end
    FigureLegend(figure[2, 1:3], history; nbanks=4, framevisible=true)
    limit = maximum(maximum(abs, filter(isfinite, vec(profile.Δκ))) for profile in differences)
    for (column, profile) in enumerate(differences)
        axis = Axis(figure[3, column]; title=profile.label * " − WRK3-SE", xlabel=L"t\;[\mathrm{days}]", ylabel=L"z\;[\mathrm{km}]")
        heatmap!(axis, profile.t / 86400, profile.z / 1000, transpose(profile.Δκ);
                 colormap=:balance, colorrange=(-limit, limit), colorscale=Symlog10(limit / 100))
    end
    Colorbar(figure[3, 4]; colormap=:balance, limits=(-limit, limit), scale=Symlog10(limit / 100), label=L"\Delta\kappa\;[\mathrm{m^2\,s^{-1}}]")
    return save_example(figure, "internal_tide_diffusivity")
end

function plot_channel(profiles, reference, differences)
    figure = Figure(size=(1050, 430))
    axes = [Axis(figure[1, column]; xlabel=label, ylabel=L"z\;[\mathrm{km}]", xscale=log10)
            for (column, label) in enumerate((L"|\kappa_h|\;[\mathrm{m^2\,s^{-1}}]", L"|\kappa|\;[\mathrm{m^2\,s^{-1}}]", L"\kappa_{\mathrm{phys}}/\kappa_{\mathrm{num}}"))]
    for profile in profiles
        color = plot_style(profile.label; distinct_colors=true).color
        for (axis, values) in zip(axes, (profile.κₕ, profile.κ, profile.ratio))
            lines!(axis, [value > 0 ? value : NaN for value in values], profile.z / 1000; color, label=profile.label)
            lines!(axis, [value < 0 ? -value : NaN for value in values], profile.z / 1000; color, linestyle=:dash)
        end
    end
    vlines!(axes[3], [1]; color=:black, linestyle=:dot)
    FigureLegend(figure[2, 1:3], axes[1]; nbanks=4, framevisible=true)
    save_example(figure, "channel_diffusivity")

    figure = Figure(size=(1250, 650))
    for (row, property, difference_property) in ((1, :Pₕ, :ΔPₕ), (2, :Pz, :ΔPz))
        values = getproperty(reference, property)
        limit = maximum(abs, filter(isfinite, vec(values)))
        axis = Axis(figure[row, 1]; title=reference.label, xlabel=L"y\;[\mathrm{km}]", ylabel=L"z\;[\mathrm{km}]")
        heatmap!(axis, reference.y / 1000, reference.z / 1000, values; colormap=:balance, colorrange=(-limit, limit), colorscale=Symlog10(limit / 100))
        Colorbar(figure[row, 2]; limits=(-limit, limit), colormap=:balance, scale=Symlog10(limit / 100))
        limit = maximum(maximum(abs, filter(isfinite, vec(getproperty(section, difference_property)))) for section in differences)
        for (column, section) in enumerate(differences)
            axis = Axis(figure[row, column+2]; title=section.label * " − WRK3-SE", xlabel=L"y\;[\mathrm{km}]", ylabel=L"z\;[\mathrm{km}]")
            heatmap!(axis, section.y / 1000, section.z / 1000, getproperty(section, difference_property);
                     colormap=:balance, colorrange=(-limit, limit), colorscale=Symlog10(limit / 100))
        end
        Colorbar(figure[row, 6]; limits=(-limit, limit), colormap=:balance, scale=Symlog10(limit / 100))
    end
    return save_example(figure, "channel_dissipation")
end

function plot_global_cost(measurements, throughput)
    figure = Figure(size=(850, 430))
    axis = Axis(figure[1, 1]; ylabel="Throughput relative to AB2-SE",
                xticks=(eachindex(measurements), [measurement.label for measurement in measurements]), xticklabelrotation=π/6)
    barplot!(axis, collect(eachindex(measurements)), throughput; color=[plot_style(measurement.label).color for measurement in measurements])
    hlines!(axis, [1]; color=:black, linestyle=:dash)
    return save_example(figure, "global_throughput")
end
