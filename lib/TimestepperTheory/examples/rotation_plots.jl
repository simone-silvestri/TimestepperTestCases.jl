using CairoMakie, LaTeXStrings

function rotation_boundary(growth, phases)
    return map(eachrow(growth)) do values
        index = findfirst(>(1e-8), values)
        isnothing(index) && return last(phases)
        index == 1 && return 0.0
        return (phases[index] + phases[index-1]) / 2
    end
end

function plot_rotation_results(results, directory)
    mkpath(directory)
    font = joinpath(pkgdir(Makie.MathTeXEngine), "assets", "fonts", "ComputerModern", "cmr10.ttf")
    set_theme!(Theme(fonts=(; regular=font), fontsize=22, linewidth=2))
    figure = Figure(size=(1380, 650), figure_padding=20)

    for (row, ε) in enumerate((0.01, 0.1)), (column, scheme) in enumerate(("WRK3", "MRK4", "SRK3"))
        names = scheme == "SRK3" ? ("SRK3-SE",) : ("$scheme-frozen-F3", "$scheme-progressive-pressure", "$scheme-SE")
        title = latexstring("\\mathrm{", scheme, "},\\quad \\varepsilon=", ε, scheme == "MRK4" ? ",\\ a=0.23" : "")
        axis = Axis(figure[row, column]; title, xlabel=L"\mu_0", ylabel=column == 1 ? L"|f|\Delta t" : "",
                    yticks=0:0.5:2, yticklabelsvisible=column == 1, yticksvisible=true)

        for (index, name) in enumerate(names)
            result = only(filter(result -> result.ε == ε && result.case.name == name, results))
            (; μ, ρ, growth, modes, θ★, case) = result
            averaging  = averaging_filter(case.averaging, case.substeps)
            limit      = (case.solver isa ForwardBackward ? 2.0 : sqrt(3.0)) / averaging.Δτ
            bounds     = rotation_boundary(growth, ρ)
            envelope   = accumulate(min, bounds)
            valid      = μ .<= min(limit, θ★ / modes.α[2])
            production = name == "$scheme-SE"

            if production
                scatter!(axis, [0.7θ★ / modes.α[2]], [scheme == "MRK4" ? 0.63 : 0.394];
                         color=:darkorange2, marker=:star5, markersize=19, strokecolor=:black, strokewidth=0.5)
                shading = band!(axis, μ[valid], zeros(count(valid)), envelope[valid]; color="#e8edf2")
                translate!(shading, 0, 0, -10)
                lines!(axis, μ[valid], bounds[valid]; color=(:black, 0.35), linewidth=1)
            end

            color = production ? :black : index == 1 ? :dodgerblue3 : :darkorange2
            style = production ? :dash  : index == 1 ? :solid       : :dashdot
            lines!(axis, μ[valid], envelope[valid]; color, linestyle=style)

            if index == 1
                θ = range(0, θ★; length=300)
                lines!(axis, θ ./ modes.α[2], sqrt.(max.(0.0, θ★^2 .- θ.^2)); color=:gray55, linestyle=:dot)
                xlims!(axis, 0, θ★ / modes.α[2])
                ylims!(axis, -0.035, 2)
                if limit < θ★ / modes.α[2]
                    vspan!(axis, limit, θ★ / modes.α[2]; color=(:gray60, 0.25))
                end
            end
        end
    end

    elements = [LineElement(color=:dodgerblue3),
                LineElement(color=:darkorange2, linestyle=:dashdot),
                LineElement(color=:black, linestyle=:dash),
                MarkerElement(color=:darkorange2, marker=:star5, markersize=19),
                LineElement(color=:gray55, linestyle=:dot)]
    labels = ["Frozen pressure and Coriolis (RK3 / F3)",
              "Progressive pressure, frozen Coriolis",
              "Progressive pressure and Coriolis (WRK3 / MRK4); frozen (SRK3)",
              "Nominal global operating point",
              "Uncoupled inertia–gravity bound"]

    Legend(figure[3, 1:3], elements, labels; nbanks=2, framevisible=true, labelsize=19)
    rowgap!(figure.layout, 12)
    colgap!(figure.layout, 16)
    save(joinpath(directory, "rotation-coupling.pdf"), figure)
    save(joinpath(directory, "rotation-coupling.png"), figure; px_per_unit=1.2)
end
