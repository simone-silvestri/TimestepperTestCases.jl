using CairoMakie
using LaTeXStrings

manuscript_font = joinpath(pkgdir(Makie.MathTeXEngine), "assets", "fonts", "ComputerModern", "cmr10.ttf")
set_theme!(Theme(fonts=(; regular=manuscript_font, bold=manuscript_font), fontsize=10, linewidth=1.6, Legend=(; framevisible=true)))

function save_results(directory, name, figure, rows)
    mkpath(directory)
    save(joinpath(directory, name * ".pdf"), figure)
    save(joinpath(directory, name * ".png"), figure; px_per_unit=220/72)
    open(joinpath(directory, name * ".csv"), "w") do stream
        println(stream, join(string.(keys(first(rows))), ','))
        for row in rows
            println(stream, join(values(row), ','))
        end
    end
    @info "Saved $name"
    return nothing
end

function plot_amplification_factor(panels, rows, output_directory)
    figure    = Figure(size=(650, 260), fontsize=11)
    transform = ReversibleScale(value -> sign(value) * log10(1 + abs(value) / 1e-10),
                                value -> sign(value) * 1e-10 * (10^abs(value) - 1))
    ticks  = [-1e-1, -1e-4, -1e-7, 0, 1e-7, 1e-4, 1e-1]
    labels = [L"-10^{-1}", L"-10^{-4}", L"-10^{-7}", L"0", L"10^{-7}", L"10^{-4}", L"10^{-1}"]
    axis   = Axis(figure[1, 1], xlabel=L"\theta", ylabel=L"|\lambda|-1", yscale=transform,
                xticks=0:0.5:3, yticks=(ticks, labels), xgridcolor=(:black, 0.1), ygridcolor=(:black, 0.1))
    for (index, panel) in enumerate(panels)
        phases, growth = panel.phases, panel.growth
        label = (L"\mathrm{AB2},\ \xi=0", L"\mathrm{AB2},\ \xi=0.1", L"\mathrm{RK2}", L"\mathrm{RK3}", L"\mathrm{RK4}", L"\mathrm{MRK4},\ a=0.23")[index]
        lines!(axis, Float64.(phases), Float64.(growth); label,
               color=(:blue, :blue, :green, :red, :purple, :purple)[index],
               linestyle=(:dash, :solid, :solid, :solid, :solid, :dash)[index])
    end
    limits!(axis, -0.05, 3.05, -1.2, 6)
    Legend(figure[1, 2], axis; framevisible=true, labelsize=10, patchsize=(24, 10), rowgap=7)
    colgap!(figure.layout, 12)
    return save_results(output_directory, "amplification_factor", figure, rows)
end

function plot_barotropic_filter_damping(panels, rows, output_directory)
    figure          = Figure(size=(600, 270), fontsize=10)
    phases          = panels[1].phases
    implicit_growth = panels[end].growth
    for (column, solver) in enumerate((ExactSolver(), ForwardBackward(), RungeKutta3()))
        axis = Axis(figure[1, column], title=("(a) Exact", "(b) FB", "(c) RK3")[column], xlabel=L"\mu_m",
                    ylabel=column == 1 ? L"|\lambda_b|-1" : "", yticklabelsvisible=column == 1,
                    xticks=0:4:16, yticks=-1:0.2:0, xgridcolor=(:black, 0.08), ygridcolor=(:black, 0.08))
        for (index, averaging) in enumerate((F1Filter(), F2Filter(), F3Filter()))
            growth = panels[3(column - 1) + index].growth
            lines!(axis, phases, growth; color=("#0072B2", "#D55E00", "#800080")[index],
                   linestyle=(:solid, :dash, :dashdot)[index], linewidth=1.3, label="F$index")
        end
        lines!(axis, phases, implicit_growth; color="#454545", linestyle=:dot, label="Implicit", linewidth=1.2)
        limits!(axis, 0, 16, -1.02, 0.025)
        column == 3 && Legend(figure[2, 1:3], axis; orientation=:horizontal, framevisible=true, patchsize=(25, 10))
    end
    colgap!(figure.layout, 10)
    return save_results(output_directory, "barotropic-filter-damping", figure, rows)
end

function plot_averaging_weights(panels, rows, output_directory)
    figure = Figure(size=(650, 350), fontsize=10)
    for (column, averaging) in enumerate((F1Filter(), F2Filter(), F3Filter()))
        (; w, Δτ, τ, w⁺, Mₚ) = panels[column]
        upper = Axis(figure[1, column], title=("F1, first order", "F2, second order", "F3, third order")[column],
                     ylabel=column == 1 ? L"w/\Delta\tau" : "", xticks=0:0.5:2, xticklabelsvisible=false,
                     yticks=([-0.5, 0, 1, 2], [L"-0.5", L"0", L"1", L"2"]),
                     yticklabelsvisible=column == 1, xgridvisible=true, xgridcolor=(:black, 0.15), xgridwidth=0.5, ygridvisible=false)
        lower = Axis(figure[2, column], xlabel=L"\tau", ylabel=column == 1 ? L"w^\dagger/\Delta\tau" : "",
                     yticks=([0, 0.5, 1], [L"0", L"0.5", L"1"]), yticklabelsvisible=column == 1,
                     xticks=([0, 0.5, 1, 1.5, 2], [L"0", L"0.5", L"1", L"1.5", L"2"]),
                     xgridvisible=true, xgridcolor=(:black, 0.15), xgridwidth=0.5, ygridvisible=false)
        barplot!(upper, τ, w ./ Δτ; width=Δτ, gap=0, color=:blue, strokewidth=0)
        barplot!(lower, τ, w⁺ ./ Δτ; width=Δτ, gap=0, color=:red, strokewidth=0)
        for (axis, limits) in ((upper, (-0.8, 2.6)), (lower, (-0.18, 1.25)))
            limits!(axis, -0.05, 2.45, limits...)
            hlines!(axis, [0]; color=:black, linewidth=0.7)
        end
        for order in 2:4
            moment = Mₚ[order - 1]
            if abs(moment) < 1e-12
                label = L"M_{%$order}=0"
            else
                exponent    = floor(Int, log10(abs(moment)))
                coefficient = round(moment / 10.0^exponent; digits=2)
                label       = L"M_{%$order}=%$coefficient\cdot10^{%$exponent}"
            end
            text!(lower, 0.97, 0.94 - 0.19 * (order - 2); text=label, space=:relative, align=(:right, :top), fontsize=9)
        end
    end
    rowgap!(figure.layout, 12)
    colgap!(figure.layout, 10)
    return save_results(output_directory, "averaging-weights", figure, rows)
end

function plot_coupling_decomposition(panels, rows, output_directory)
    figure = Figure(size=(576, 331), fontsize=12)
    for (row, ε) in enumerate((big"0.0005", big"0.05")), (column, scheme) in enumerate((WickerSkamarock3(), ShuOsher3()))
        (; phases, corrections, complete, retained, predicted) = panels[2(row - 1) + column]
        axis = Axis(figure[row, column], xlabel=L"\mu_0", title=row == 1 ? (scheme isa ShuOsher3 ? "SRK3" : "WRK3") : "")
        lines!(axis, Float64.(phases), Float64.([result.growth / result.damping for result in corrections]); color=:black, label="Full matrix")
        lines!(axis, Float64.(phases), Float64.(complete); color="#CC79A7", linestyle=:dash, label="All second-order terms")
        scatter!(axis, Float64.(phases[1:5:end]), Float64.(retained[1:5:end]); color=:white, strokecolor="#0072B2",
                 strokewidth=0.8, markersize=4, label="Direct + barotropic indirect")
        lines!(axis, Float64.(phases), Float64.(predicted); color="#D55E00", linestyle=:dashdot, label="Leading analytical prediction")
        text!(axis, 0.03, 0.96; space=:relative, align=(:left, :top), text=row == 1 ? L"\varepsilon=5\times10^{-4}" : L"\varepsilon=0.05")
        xlims!(axis, 0, 16)
        row == 2 && column == 2 && Legend(figure[3, 1:2], axis; orientation=:horizontal, nbanks=2, labelsize=10.5, patchsize=(20, 10), rowgap=3, padding=(5, 5, 5, 5))
    end
    Label(figure[1:2, 0], L"\mathrm{Growth}/(1-|\lambda_j|)"; rotation=π/2)
    return save_results(output_directory, "coupling_decomposition", figure, rows)
end

function plot_damping_balance(panels, rows, output_directory)
    figure = Figure(size=(626, 468), fontsize=12)
    for (column, scheme) in enumerate((WickerSkamarock3(), ShuOsher3()))
        axes = [Axis(figure[row, column], xlabel=L"\mu_0", title=row == 1 ? (scheme isa ShuOsher3 ? "SRK3" : "WRK3") : "",
                     ylabel=column == 1 ? ("Direct", "Barotropic indirect", "Total growth")[row] : "") for row in 1:3]
        configurations = ((ExactSolver(), NoFilter()), (ExactSolver(), F3Filter()), (RungeKutta3(), F3Filter()),
                          (ForwardBackward(), F3Filter()), (ForwardBackward(), F1Filter()))
        for (index, (solver, averaging)) in enumerate(configurations)
            (; phases, measured, predicted) = panels[5(column - 1) + index]
            color     = ("#666666", "#AA4499", "#0072B2", "#D55E00", "#009E73")[index]
            label     = ("Exact, unfiltered", "Exact + F3", "RK3 + F3", "FB + F3", "FB + F1")[index]
            linestyle = (:solid, :dash, :solid, :dash, :dashdot)[index]
            for row in 1:3
                lines!(axes[row], Float64.(phases), Float64.(measured[row]); label, color=(color, 0.72), linestyle, linewidth=1.1)
                targets = [exp10.(range(log10(0.05), 0; length=9)) .* (1 + 0.035 * (index - 1));
                           sort(unique([collect(1.4:0.65:20); collect(4.8:0.3:7.1)])) .+ 0.07 * (index - 1)]
                indices = unique([argmin(abs.(phases .- target)) for target in targets])
                scatter!(axes[row], Float64.(phases[indices]), Float64.(predicted[row][indices]); color=:white,
                         strokecolor=color, strokewidth=0.7, markersize=3.4)
            end
        end
        for row in 1:3
            limits = ((-1.08, 0.22), (-0.55, 0.7), (-2.5, 0.5))[row]
            row == 1 && (axes[row].yticks = -1:0.25:0)
            row == 2 && (axes[row].yticks = -0.5:0.25:0.5)
            limits!(axes[row], 0, 20, limits...)
            hlines!(axes[row], [0]; color=:gray, linewidth=0.6)
            row == 3 && hlines!(axes[row], [-1]; color=:gray, linestyle=:dot, linewidth=0.7)
        end
        if column == 2
            lines!(axes[3], [NaN], [NaN]; color=(:gray, 0.72), label=L"\mathrm{Matrix},\ \varepsilon=0.01")
            scatter!(axes[3], [NaN], [NaN]; color=:white, strokecolor=:black, strokewidth=0.7, markersize=4, label="Leading approximation")
            Legend(figure[4, 1:2], axes[3]; orientation=:horizontal, nbanks=3, labelsize=10.5, patchsize=(20, 8), rowgap=4, padding=(5, 5, 5, 5))
        end
    end
    Label(figure[1:3, 0], L"\mathrm{Contribution}/(1-|\lambda_j|)"; rotation=π/2)
    return save_results(output_directory, "damping_balance", figure, rows)
end

function plot_forcing_stratification(panels, rows, output_directory)
    figure = Figure(size=(950, 460), fontsize=15)
    ticks  = (10.0 .^ (-6:2:0), [L"10^{-6}", L"10^{-4}", L"10^{-2}", L"1"])
    axes   = [Axis(figure[1, column], xscale=log10, yscale=log10, xticks=ticks, xticksvisible=false, xticklabelsvisible=false,
                 ylabel=column == 1 ? "Absolute magnitude\ncorrection" : "", title=("Direct", "Barotropic indirect", "Combined correction")[column]) for column in 1:3]
    ratios = [Axis(figure[2, column], xlabel=column == 2 ? L"\varepsilon" : "", ylabel=column == 1 ? "Computed /\npredicted" : "",
                   xscale=log10, yticks=WilkinsonTicks(3), xticks=ticks) for column in 1:3]
    linkyaxes!(axes...)
    for column in 1:3
        linkxaxes!(axes[column], ratios[column])
        xlims!(axes[column], 4e-7, 1)
        hlines!(ratios[column], [1]; color=:gray, linestyle=:dot)
    end
    for (index, panel) in enumerate(panels)
        (; forcing, stratifications, measured, predicted) = panel
        color = ("#0072B2", "#009E73")[index]
        label = ("Frozen", "Progressive")[index]
        for column in 1:3
            lines!(axes[column], Float64.(stratifications), Float64.(abs.(measured[column])); label, color, linewidth=1.4)
            if forcing isa FrozenForcing
                scatter!(axes[column], Float64.(stratifications), Float64.(abs.(predicted[column])); marker=:circle,
                         color=:white, strokecolor=color, strokewidth=1.1, markersize=6)
                lines!(ratios[column], Float64.(stratifications), Float64.(measured[column] ./ predicted[column]); color, linewidth=1.4)
            end
            exponent  = forcing isa FrozenForcing ? 2 : 3
            guide     = reverse(stratifications)[3:5]
            amplitude = abs(reverse(measured[column])[4]) * 8 / guide[2]^exponent
            lines!(axes[column], Float64.(guide), Float64.(amplitude .* guide.^exponent); color=:gray, linestyle=:dot)
            text!(axes[column], Float64(guide[end]), Float64(amplitude * guide[end]^exponent);
                  text=exponent == 2 ? L"\varepsilon^2" : L"\varepsilon^3", align=(:left, :bottom))
        end
    end
    scatter!(axes[1], [NaN], [NaN]; marker=:circle, color=:white, strokecolor=:black, strokewidth=1.1,
             markersize=6, label="Frozen prediction")
    rowgap!(figure.layout, 1, 12)
    rowsize!(figure.layout, 2, Relative(0.23))
    Legend(figure[3, 1:3], axes[1]; orientation=:horizontal, labelsize=14, patchsize=(26, 10),
           rowgap=2, padding=(7, 7, 7, 7))
    return save_results(output_directory, "forcing_stratification", figure, rows)
end

function plot_mrk4_forcing(panels, rows, output_directory)
    figure = Figure(size=(1100, 500), fontsize=20)
    axes   = Axis[]
    for (column, panel) in enumerate(panels)
        (; ε, phases, results, prediction_phases, predictions, damping) = panel
        title = column == 1 ? L"\varepsilon=5\cdot10^{-4},\quad a=1/4" : L"\varepsilon=0.05,\quad a=1/4"
        axis  = Axis(figure[1, column]; title, xlabel=L"\mu_0", xscale=Makie.Symlog10(5; linscale=0.6),
                    xticks=[0, 2, 5, 10, 20, 30], xgridcolor=(:black, 0.12), ygridcolor=(:black, 0.12))
        column == 1 && (axis.ylabel = L"\mathrm{Coupling\ correction}\,/\,\theta_1^4")
        push!(axes, axis)
        for (index, label, color) in ((1, "Frozen", "#0072B2"), (2, "Progressive", "#009E73"))
            lines!(axis, Float64.(phases), getindex.(results, index); label, color, linewidth=2)
        end
        scatter!(axis, Float64.(prediction_phases), predictions; label="Frozen prediction",
                 marker=:circle, color=:transparent, strokecolor="#0072B2", strokewidth=1.2, markersize=7)
        lines!(axis, Float64.(phases), getindex.(results, 3); color=:black, linestyle=:dash, linewidth=2, label="Damping (a = 1/4)")
        for (index, color, linestyle, label) in ((1, :black, :dot, "Damping (a = 0)"), (2, :darkorange3, :dashdot, "Damping (a = 0.23)"))
            lines!(axis, Float64.(phases), Float64.(damping[index]); color, linestyle, label, linewidth=2)
        end
        xlims!(axis, 0, 30)

        inset = Axis(figure[1, column], width=Relative(0.29), height=Relative(0.23), halign=0.16, valign=0.72,
                     title="Small-step limit", titlesize=14, xlabel=L"\mu_0", xlabelsize=17, xticklabelsize=13, yticklabelsize=13,
                     xticks=[0, 0.5, 1], backgroundcolor=:white, xgridcolor=(:black, 0.08), ygridcolor=(:black, 0.08))
        translate!(inset.blockscene, 0, 0, 100)
        selected = findall(<=(1), phases)
        lines!(inset, Float64.(phases[selected]), getindex.(results[selected], 2); color="#009E73", linewidth=1.5)
        lines!(inset, Float64.(phases[selected]), getindex.(results[selected], 3); color=:black, linestyle=:dash, linewidth=1.2)
        hlines!(inset, [0]; color=(:gray, 0.5), linewidth=0.6)
        xlims!(inset, -0.05, 1.05)
        ylims!(inset, column == 1 ? (-7e-5, 1.5e-5) : (-1.5e-5, 1e-4))
        inset.yticks = column == 1 ? ([-6e-5, -3e-5, 0], [L"-6\cdot10^{-5}", L"-3\cdot10^{-5}", L"0"]) :
                                    ([0, 4e-5, 8e-5], [L"0", L"4\cdot10^{-5}", L"8\cdot10^{-5}"])
    end
    linkyaxes!(axes...)
    Legend(figure[2, 1:2], axes[1]; orientation=:horizontal, nbanks=2, labelsize=17,
           rowgap=3, colgap=16, patchsize=(28, 10), padding=(8, 8, 5, 5))
    return save_results(output_directory, "mrk4_forcing", figure, rows)
end

function plot_schemes(panels, rows, output_directory)
    figure = Figure(size=(540, 390), fontsize=10, figure_padding=(12, 16, 10, 14))
    names  = ["AB2-SE", "WRK3-SE-F1", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM"]
    for (column, count) in enumerate((8, 32)), (row, ε) in enumerate((big"0.05", big"0.15"))
        (; phases, growth, thresholds) = panels[2(column - 1) + row]
        axis = Axis(figure[row, column], xlabel=row == 2 ? L"\mu_0" : "", ylabel=column == 1 ? L"\max |\lambda| - 1" : "",
                    xticklabelsvisible=row == 2, yticklabelsvisible=column == 1,
                    xticks=0:10:40, titlegap=8, yticklabelspace=39.0,
                    xgridcolor=(:black, 0.08), ygridcolor=(:black, 0.08),
                    yscale=Makie.Symlog10(1e-13; linscale=0.35 * 10/9),
                    yticks=([-1e-2, -1e-6, -1e-10, 0, 1e-6], [L"-10^{-2}", L"-10^{-6}", L"-10^{-10}", L"0", L"10^{-6}"]),
                    title=L"\varepsilon=%$(Float64(ε)),\ N_m=%$count")
        limits!(axis, 0, 44, -6e-2, 3e-4)
        for (index, name) in enumerate(names)
            lines!(axis, Float64.(phases), growth[index]; label=name,
                   color=("#0000ff", "#ff0000", "#ff0000", "#008000", "#800080", "#ff0000")[index],
                   linestyle=(:solid, :dot, :solid, :dash, :solid, :dashdot)[index], linewidth=index == 2 ? 1.1 : 0.9)
            @info "Stability" scheme=name modes=count stratification=Float64(ε)
        end
        vlines!(axis, thresholds; color=:gray, linestyle=:dot)
        row == 2 && column == 2 && Legend(figure[3, 1:2], axis; nbanks=2, orientation=:horizontal, labelsize=9, patchsize=(25, 10), rowgap=5, colgap=18, padding=(8, 8, 6, 6))
    end
    rowgap!(figure.layout, 1, 6)
    rowgap!(figure.layout, 2, 14)
    colgap!(figure.layout, 12)
    return save_results(output_directory, "schemes", figure, rows)
end

function plot_order_schemes(panels, rows, output_directory)
    figure = Figure(size=(461, 220), fontsize=9)
    axis   = Axis(figure; bbox=BBox(50.7, 336.5, 39.6, 213.4), xlabel=L"\mu_0", ylabel="Global error", xscale=log10, yscale=log10, xlabelsize=10, ylabelsize=10)
    limits!(axis, 8e-4, 0.6, 1e-10, 4)
    axis.yticks = (10.0 .^ (-10:2:0), [L"10^{%$exponent}" for exponent in -10:2:0])
    axis.xgridvisible = false
    for (index, name) in enumerate(["AB2-SE", "WRK3-SE-F1", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM"])
        results = panels[index]
        color   = ("#0000ff", "#ff0000", "#ff0000", "#008000", "#800080", "#ff0000")[index]
        lines!(axis, Float64.([result.phase for result in results]), [result.error for result in results]; label=name, color,
               linestyle=(:solid, :dot, :solid, :dash, :solid, :dashdot)[index])
        if index in (2, 5)
            scatter!(axis, Float64.([result.phase for result in results][2:4:end]), [result.error for result in results][2:4:end];
                     marker=index == 2 ? :rect : :diamond, color=:white, strokecolor=color, strokewidth=0.9, markersize=4, label=name)
        end
    end
    for exponent in 1:2
        guide     = [exponent == 1 ? 0.0015 : 0.004, 0.025]
        reference = (exponent == 1 ? 0.18 : 1e-5) .* guide.^exponent
        lines!(axis, guide, reference; color=:gray, linestyle=:dash, linewidth=0.8)
        text!(axis, guide[end], reference[end]; offset=(3, 4), align=(:left, :bottom), fontsize=8, color=:gray,
              text=exponent == 1 ? L"\mathrm{First\ order}:\ \mu_0" : L"\mathrm{Second\ order}:\ \mu_0^2")
    end
    Legend(figure, axis; bbox=BBox(348, 461, 45, 205), merge=true, labelsize=9, rowgap=6, patchsize=(20, 9), padding=(7, 7, 7, 7))
    return save_results(output_directory, "order_schemes", figure, rows)
end
