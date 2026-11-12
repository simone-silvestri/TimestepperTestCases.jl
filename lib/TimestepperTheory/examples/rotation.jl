# # Coriolis forcing in split-explicit stages
# This example reproduces Supporting Information S8 using the analytical modal timestep matrices.
# The stationary geostrophic mode contributes unit amplification; plotted positive values indicate instability.

using TimestepperTheory
using TimestepperTheory: rotation_outputs, stage_fractions
using LinearAlgebra
using CairoMakie
using LaTeXStrings

function write_rotation_results(directory)
    couplings = (FrozenRotation(), ProgressivePressure(), ProgressiveRotation(), SSPRotation(), SSPSubstepRotation())
    labels = ("WRK3: frozen", "WRK3: progressive pressure, frozen Coriolis", "WRK3: progressive pressure and Coriolis",
              "SRK3: frozen pressure and Coriolis", "SRK3: substepped Coriolis")
    colors = (:black, :dodgerblue3, :seagreen3, :darkorange2, :purple)
    styles = (:solid, :dash, :dot, :dashdot, :solid)
    μ = collect(range(0.05, 6., length=240))
    ρ = collect(range(0., 2., length=161))
    filter = averaging_filter(F3Filter(), 96)
    font = joinpath(pkgdir(Makie.MathTeXEngine), "assets", "fonts", "ComputerModern", "cmr10.ttf")
    set_theme!(Theme(fonts=(; regular=font, bold=font), fontsize=18, linewidth=2.3))
    figure = Figure(size=(950, 420), figure_padding=14)
    transform = ReversibleScale(x -> sign(x) * log10(1 + abs(x) / 1e-10), x -> sign(x) * 1e-10 * (10^abs(x) - 1))
    ticks = ([-1e-7, 0, 1e-7, 1e-4, 1.],
             [L"-10^{-7}", L"0", L"10^{-7}", L"10^{-4}", L"1"])
    rows = NamedTuple[]
    thresholds = NamedTuple[]
    for (column, ε) in enumerate((0.01, 0.1))
        modes = vertical_modes(ε, 8)
        bottom = Axis(figure[1, column], title=latexstring("\\varepsilon = $ε"), xlabel=L"f\Delta t",
                      ylabel=column == 1 ? L"\max_{\mu_0,j}|\lambda_j|-1" : "", yscale=transform, yticks=ticks)
        for (index, coupling) in enumerate(couplings)
            scheme = WickerSkamarock3()
            filter = averaging_filter(F3Filter(), 96)
            c = coupling isa Union{SSPRotation, SSPSubstepRotation} ? (1., 1., 1.) : stage_fractions(scheme, Float64)
            outputs = [[rotation_outputs(μᵢ, cᵢ, filter) for cᵢ in c] for μᵢ in μ]
            growth = zeros(length(ρ))
            locations = similar(growth)
            for (i, ρᵢ) in enumerate(ρ)
                values = [maximum(abs, eigvals(rotating_evolution_matrix(μᵢ, ρᵢ, modes, coupling; scheme, filter, outputs=outputs[j]))) - 1 for (j, μᵢ) in enumerate(μ)]
                growth[i], j = findmax(values)
                locations[i] = μ[j]
                push!(rows, (; ε, scheme=labels[index], ρ=ρᵢ, μ=μ[j], growth=growth[i]))
            end
            lines!(bottom, ρ, growth; color=colors[index], linestyle=styles[index], label=labels[index])
            first = findfirst(>(1e-8), growth)
            localgrowth = maximum(abs, eigvals(rotating_evolution_matrix(1.58, 0.5, modes, coupling; scheme, filter))) - 1
            push!(thresholds, (; ε, scheme=labels[index], growth=localgrowth,
                               threshold=isnothing(first) ? NaN : ρ[first], μ=isnothing(first) ? NaN : locations[first]))
            @info "Completed" ε scheme=labels[index] threshold=last(thresholds)
        end
        xlims!(bottom, 0, 2)
        ylims!(bottom, -1e-7, 4)
    end
    Legend(figure[2, 1:2], figure.content[1]; nbanks=2, framevisible=true, labelsize=15, patchsize=(30, 12), rowgap=5, padding=(8, 8, 6, 6))
    rowgap!(figure.layout, 10)
    save(joinpath(directory, "rotation-coupling.pdf"), figure)
    save(joinpath(directory, "rotation-coupling.png"), figure; px_per_unit=2)
    for (name, data) in (("rotation_sweep.csv", rows), ("rotation_table.csv", thresholds))
        open(joinpath(directory, name), "w") do io
            println(io, join(keys(first(data)), ','))
            for row in data
                println(io, join(values(row), ','))
            end
        end
    end
    return thresholds
end

if abspath(PROGRAM_FILE) == @__FILE__
    write_rotation_results(mkpath(joinpath(@__DIR__, "..", "output", "rotation")))
end
