# # Time discretization of free-surface ocean models
#
# This example reproduces the amplification, filtering, and coupled modal analyses of Sections 2–4. The numerical implementation is provided by
# TimestepperTheory. Each experiment shows its computations; figure formatting and file output are hidden in the rendered page.

using TimestepperTheory
using CairoMakie #src
using LaTeXStrings #src
using LinearAlgebra

include("time_discretization_figures.jl") #src

manuscript_font = joinpath(pkgdir(Makie.MathTeXEngine), "assets", "fonts", "ComputerModern", "cmr10.ttf") #src
set_theme!(Theme(fonts=(; regular=manuscript_font, bold=manuscript_font), fontsize=10, linewidth=1.6, #src
                 Legend=(; framevisible=true))) #src

output_directory = joinpath(@__DIR__, "..", "output")
samples = "--quick" in ARGS ? 9 : 61
setprecision(256)

# ## Amplification factors of the outer schemes
#
# For a Fourier mode of unit-speed advection, the exact timestep multiplies the state by exp(−iθ), where θ = kΔt.
# The magnitude departure |λ| − 1 measures numerical damping or growth. WRK3 and SRK3 share the same cubic polynomial;
# RK2 and the fourth-order Taylor polynomial provide references. The AB2 curves use the physical root, with ξ = 0 and 0.1.

let
    figure = Figure(size=(650, 260), fontsize=11) #src
    transform = ReversibleScale(value -> sign(value) * log10(1 + abs(value) / 1e-10), #src
                                value -> sign(value) * 1e-10 * (10^abs(value) - 1)) #src
    ticks = [-1e-1, -1e-4, -1e-7, 0, 1e-7, 1e-4, 1e-1] #src
    labels = [L"-10^{-1}", L"-10^{-4}", L"-10^{-7}", L"0", L"10^{-7}", L"10^{-4}", L"10^{-1}"] #src
    axis = Axis(figure[1, 1], xlabel=L"\theta", ylabel=L"|\lambda|-1", yscale=transform, #src
                xticks=0:0.5:3, yticks=(ticks, labels), xgridcolor=(:black, 0.1), ygridcolor=(:black, 0.1)) #src
    rows = NamedTuple[]
    phases = BigFloat.(0:300) ./ 100
    for (index, ξ) in enumerate((0, 1 // 10))
        growth = [abs(amplification_factor(phase, AdamsBashforth2(); ξ)) - 1 for phase in phases]
        label = index == 1 ? L"\mathrm{AB2},\ \xi=0" : L"\mathrm{AB2},\ \xi=0.1" #src
        lines!(axis, Float64.(phases), Float64.(growth); color=:blue, linestyle=index == 1 ? :dash : :solid, label) #src
        append!(rows, [(; scheme="AB2", parameter=Float64(ξ), phase, growth=value) for (phase, value) in zip(phases, growth)])
    end
    for (index, scheme) in enumerate((RungeKutta2(), WickerSkamarock3(), ModifiedRungeKutta4(a=1 // 4), ModifiedRungeKutta4()))
        growth = [abs(amplification_factor(phase, scheme)) - 1 for phase in phases]
        label = (L"\mathrm{RK2}", L"\mathrm{RK3}", L"\mathrm{RK4}", L"\mathrm{MRK4},\ a=0.23")[index] #src
        lines!(axis, Float64.(phases), Float64.(growth); color=(:green, :red, :purple, :purple)[index], #src
               linestyle=index == 4 ? :dash : :solid, label) #src
        name = ("RK2", "RK3", "RK4", "MRK4")[index]
        append!(rows, [(; scheme=name, parameter=scheme isa ModifiedRungeKutta4 ? Float64(scheme.a) : NaN, phase, growth=value) for (phase, value) in zip(phases, growth)])
    end
    limits!(axis, -0.05, 3.05, -1.2, 6) #src
    Legend(figure[1, 2], axis; framevisible=true, labelsize=10, patchsize=(24, 10), rowgap=7) #src
    colgap!(figure.layout, 12) #src
    save_results(output_directory, "amplification_factor", figure, rows) #src
end

#md # ![](amplification_factor.png)

# ## Barotropic substeps and averaging filters
#
# A barotropic solve averages the substep evolution with weights w. The three panels compare exact, forward–backward,
# and RK3 substeps, using the F1, F2, and F3 filters and 64 nominal substeps. The grey curve is the implicit solver.

let
    figure = Figure(size=(600, 270), fontsize=10) #src
    phases = sort(unique([0; exp10.(range(-2, 0; length=180)); range(1, 16; length=901)]))
    implicit_growth = @. inv(sqrt(1 + phases^2)) - 1
    rows = NamedTuple[]
    for (column, solver) in enumerate((ExactSolver(), ForwardBackward(), RungeKutta3()))
        axis = Axis(figure[1, column], title=("(a) Exact", "(b) FB", "(c) RK3")[column], xlabel=L"\mu_m", #src
                    ylabel=column == 1 ? L"|\lambda_b|-1" : "", yticklabelsvisible=column == 1, #src
                    xticks=0:4:16, yticks=-1:0.2:0, xgridcolor=(:black, 0.08), ygridcolor=(:black, 0.08)) #src
        for (index, averaging) in enumerate((F1Filter(), F2Filter(), F3Filter()))
            filter = averaging_filter(averaging, 64)
            growth = map(phases) do phase
                response = barotropic_response(phase; solver, filter)
                return sum(abs, eigvals(response.W)) / 2 - 1
            end
            lines!(axis, phases, growth; color=("#0072B2", "#D55E00", "#800080")[index], #src
                   linestyle=(:solid, :dash, :dashdot)[index], linewidth=1.3, label="F$index") #src
            name = ("Exact", "FB", "RK3")[column]
            append!(rows, [(; solver=name, filter="F$index", phase, growth=value) for (phase, value) in zip(phases, growth)])
        end
        lines!(axis, phases, implicit_growth; color="#454545", linestyle=:dot, label="Implicit", linewidth=1.2) #src
        limits!(axis, 0, 16, -1.02, 0.025) #src
        column == 3 && Legend(figure[2, 1:3], axis; orientation=:horizontal, framevisible=true, patchsize=(25, 10)) #src
    end
    append!(rows, [(; solver="Implicit", filter="none", phase, growth=value) for (phase, value) in zip(phases, implicit_growth)])
    colgap!(figure.layout, 10) #src
    save_results(output_directory, "barotropic-filter-damping", figure, rows) #src
end

#md # ![](barotropic-filter-damping.png)

# The transport weights follow by summing the discrete continuity updates: each weight is the fractional spacing
# multiplied by the sum of all subsequent output weights, including the current one. Both sets of weights integrate
# to one. The centered moments M₂, M₃, and M₄ describe the filter's low-frequency damping and phase error.

let
    figure = Figure(size=(650, 350), fontsize=10) #src
    rows = NamedTuple[]
    for (column, averaging) in enumerate((F1Filter(), F2Filter(), F3Filter()))
        filter = averaging_filter(averaging, 64)
        times = eachindex(filter.w) .* filter.Δτ
        transport = filter.Δτ .* reverse(cumsum(reverse(filter.w)))
        upper = Axis(figure[1, column], title=("F1, first order", "F2, second order", "F3, third order")[column], #src
                     ylabel=column == 1 ? L"w/\Delta\tau" : "", xticks=0:0.5:2, xticklabelsvisible=false, #src
                     yticks=([-0.5, 0, 1, 2], [L"-0.5", L"0", L"1", L"2"]), #src
                     yticklabelsvisible=column == 1, xgridvisible=false, ygridvisible=false) #src
        lower = Axis(figure[2, column], xlabel=L"\tau", ylabel=column == 1 ? L"w^\dagger/\Delta\tau" : "", #src
                     yticks=([0, 0.5, 1], [L"0", L"0.5", L"1"]), yticklabelsvisible=column == 1, #src
                     xticks=([0, 0.5, 1, 1.5, 2], [L"0", L"0.5", L"1", L"1.5", L"2"]), #src
                     xgridvisible=false, ygridvisible=false) #src
        barplot!(upper, times, filter.w ./ filter.Δτ; width=filter.Δτ, gap=0, color=:blue, strokewidth=0) #src
        barplot!(lower, times, transport ./ filter.Δτ; width=filter.Δτ, gap=0, color=:red, strokewidth=0) #src
        for (axis, limits) in ((upper, (-0.8, 2.6)), (lower, (-0.18, 1.25))) #src
            limits!(axis, -0.05, 2.45, limits...) #src
            hlines!(axis, [0]; color=:black, linewidth=0.7) #src
            vlines!(axis, [0, 1]; color=(:black, 0.45), linewidth=0.5) #src
        end #src
        for order in 2:4
            moment = sum(filter.w .* (times .- 1).^order)
            if abs(moment) < 1e-12 #src
                label = L"M_{%$order}=0" #src
            else #src
                exponent = floor(Int, log10(abs(moment))) #src
                coefficient = round(moment / 10.0^exponent; digits=2) #src
                label = L"M_{%$order}=%$coefficient\cdot10^{%$exponent}" #src
            end #src
            text!(lower, 0.97, 0.94 - 0.19 * (order - 2); text=label, space=:relative, align=(:right, :top), fontsize=9) #src
        end
        append!(rows, [(; filter="F$column", time, spacing=filter.Δτ, weight, transport_weight)
                       for (time, weight, transport_weight) in zip(times, filter.w, transport)])
    end
    rowgap!(figure.layout, 12) #src
    colgap!(figure.layout, 10) #src
    save_results(output_directory, "averaging-weights", figure, rows) #src
end

#md # ![](averaging-weights.png)

# ## Vertical modes and the modal state
#
# The dimensionless stratification is ε = N²H/g. Each root of x tan(x) = ε determines a normalized wave speed α,
# a surface value ψ, and a depth average χ. The requested count excludes the barotropic mode, which is always retained.

stratification = big"0.05"
modes = vertical_modes(stratification, 2)

modes.α

# The state is ordered as y = [u₀, …, uQ, h₀, …, hQ], with velocity divided by √(gH) and pressure divided by gH.
# The depth average is obtained by multiplying the corresponding amplitude vector by χᵀ. At finite mode count,
# χᵀχ is slightly less than one. The modal weights retain this truncation error; the SRK3 blended matrices use χᵀχ ≈ 1.

omitted_weight = 1 - sum(abs2, modes.χ)

# ## Barotropic integration and polynomial forcing
#
# The barotropic Courant number μ₀ = k√(gH)Δt uses the reference speed √(gH). The averaging filter
# determines the returned velocity and continuity transport. F3 has vanishing second and third centered moments.

phase = big"2.0"
filter = averaging_filter(F3Filter(); float_type=BigFloat)
response = barotropic_response(phase; solver=RungeKutta3(), filter)

response.W

# The forcing outputs correspond to 1, s, and s², where s is time divided by the stage duration. Multiplying these
# outputs by the reconstructed polynomial coefficients gives the forced velocity and continuity transport.

response.ūₚ, response.ūₚ⁺

# The continuity transport follows from the same discrete surface update, including the averaging filter.

continuity_residual = norm(response.W[2, :] - [0, 1] + im * phase * response.W̃)

# ## Stage matrices and the complete timestep
#
# For WRK3-SE, the final stage uses a quadratic reconstruction of the slow forcing. The returned arrays contain
# A[m] = Aₘ, B[m] = Bₘ, and M[m] = Mₘ. The initial matrix M₀ = I is used internally. Contributions from earlier stages enter through Aₘ.

stages = stage_matrices(phase, modes; scheme=WickerSkamarock3(), forcing=FinalStageQuadraticForcing(), solver=RungeKutta3(), filter)
Aₘ = stages.A[end]
Bₘ = stages.B[end]
Mₘ = stages.M[end]

Mₘ

# SRK3 uses its blended stage update and SSP-weighted forcing. MRK4 uses progressive reconstruction, including a
# least-squares quadratic at the fourth stage. The implicit WRK3 scheme uses an implicit surface solve. AB2 uses an
# augmented state [yⁿ; yⁿ⁻¹], so its complete matrix has twice the dimension of a Runge–Kutta timestep matrix.

configurations = ("AB2-SE", "WRK3-SE-F1", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM")

for name in configurations
    configuration = scheme_configuration(name)
    averaging = averaging_filter(configuration.averaging; float_type=BigFloat)
    matrix = evolution_matrix(phase, modes; configuration.scheme, configuration.forcing, configuration.solver, filter=averaging)
    println(name, ": matrix size ", size(matrix), ", maximum growth ", Float64(maximum(abs, eigvals(matrix)) - 1))
end

# ## Direct and indirect coupling corrections
#
# The uncoupled matrix advances each mode independently. Projecting the difference between the complete and uncoupled
# matrices onto their normalized left and right eigenvectors gives the coupling coefficients Vⱼₗ. These determine the
# direct correction within the selected mode and propagation direction, and the indirect correction through other modes.

reference_matrix = uncoupled_matrix(phase, modes; scheme=WickerSkamarock3(), solver=RungeKutta3(), filter)
corrections = coupling_corrections(complete_matrix, reference_matrix, modes)

corrections.direct_magnitude, corrections.indirect_magnitude, corrections.other_magnitude

# The first two figures compare the perturbation expansion with the full eigenvalue and separate the damping balance
# for exact, RK3, and forward–backward barotropic integration. These calculations use 128-bit arithmetic and eight
# baroclinic modes. The decomposition includes both stratifications shown in the manuscript figure.

setprecision(128) do
    filter = averaging_filter(F3Filter(); float_type=BigFloat)
    phases = samples == 9 ? BigFloat.(range(0.1, 16; length=9)) : BigFloat.(1:160) ./ 10
    figure = Figure(size=(576, 331), fontsize=10) #src
    rows = NamedTuple[]
    for (row, stratification) in enumerate((big"0.0005", big"0.05")), (column, scheme) in enumerate((WickerSkamarock3(), ShuOsher3()))
        forcing = scheme == ShuOsher3() ? SSPForcing() : FrozenForcing()
        axis = Axis(figure[row, column], xlabel=L"\mu_0", title=row == 1 ? (scheme isa ShuOsher3 ? "SRK3" : "WRK3") : "") #src
        modes = vertical_modes(stratification, 8)
        corrections = [coupling_corrections(evolution_matrix(phase, modes; scheme, forcing, filter),
                                            uncoupled_matrix(phase, modes; scheme, filter), modes) for phase in phases]
        predictions = [frozen_predictions(phase, modes; scheme, filter) for phase in phases]
        complete = [-1 + (result.direct_magnitude + result.indirect_magnitude + result.other_magnitude) / result.damping for result in corrections]
        retained = [-1 + (result.direct_magnitude + result.indirect_magnitude) / result.damping for result in corrections]
        predicted = [-1 + real(result.reference_eigenvalue) / abs(result.reference_eigenvalue) * prediction.total / result.damping for (prediction, result) in zip(predictions, corrections)]
        lines!(axis, Float64.(phases), Float64.([result.growth / result.damping for result in corrections]); color=:black, label="Full matrix") #src
        lines!(axis, Float64.(phases), Float64.(complete); color="#CC79A7", linestyle=:dash, label="All second-order terms") #src
        scatter!(axis, Float64.(phases[1:5:end]), Float64.(retained[1:5:end]); color=:white, strokecolor="#0072B2", #src
                 strokewidth=0.8, markersize=4, label="Direct + barotropic indirect") #src
        lines!(axis, Float64.(phases), Float64.(predicted); color="#D55E00", linestyle=:dashdot, label="Leading analytical prediction") #src
        text!(axis, 0.03, 0.96; space=:relative, align=(:left, :top), text=row == 1 ? L"\varepsilon=5\times10^{-4}" : L"\varepsilon=0.05") #src
        xlims!(axis, 0, 16) #src
        for (index, phase) in enumerate(phases)
            result = corrections[index]
            push!(rows, (; scheme, phase, stratification, direct=result.direct_magnitude, indirect=result.indirect_magnitude,
                          other=result.other_magnitude, growth=result.growth, damping=result.damping, prediction=predicted[index], overlap=result.overlap))
        end
        row == 2 && column == 2 && Legend(figure[3, 1:2], axis; orientation=:horizontal, nbanks=2, labelsize=8.5, patchsize=(20, 10), rowgap=3, padding=(5, 5, 5, 5)) #src
    end
    Label(figure[1:2, 0], L"\mathrm{Growth}/(1-|\lambda_j|)"; rotation=π/2) #src
    save_results(output_directory, "coupling_decomposition", figure, rows) #src
end

#md # ![](coupling_decomposition.png)

# ### Damping balance
#
# We repeat the decomposition with each barotropic solver and filter, and compare its direct and indirect contributions.

setprecision(128) do
    figure = Figure(size=(626, 468), fontsize=9) #src
    rows = NamedTuple[]
    for (column, scheme) in enumerate((WickerSkamarock3(), ShuOsher3()))
        forcing = scheme == ShuOsher3() ? SSPForcing() : FrozenForcing()
        axes = [Axis(figure[row, column], xlabel=L"\mu_0", title=row == 1 ? (scheme isa ShuOsher3 ? "SRK3" : "WRK3") : "", #src
                     ylabel=column == 1 ? ("Direct", "Barotropic indirect", "Total growth")[row] : "") for row in 1:3] #src
        configurations = ((ExactSolver(), NoFilter()), (ExactSolver(), F3Filter()), (RungeKutta3(), F3Filter()),
                          (ForwardBackward(), F3Filter()), (ForwardBackward(), F1Filter()))
        for (index, (solver, averaging)) in enumerate(configurations)
            response_name = ("exact", "exact_F3", "RK3_F3", "FB_F3", "FB_F1")[index]
            averaging = averaging_filter(averaging; float_type=BigFloat)
            modes = vertical_modes(big"0.01", 8)
            phases = samples == 9 ? BigFloat.(range(0.05, 20; length=9)) : BigFloat.(sort(unique([collect(0.05:0.025:1); collect(1:0.1:20)])))
            if response_name == "exact" && samples != 9
                append!(phases, [2BigFloat(π)*resonance + BigFloat(offset)/500 for resonance in 1:3 for offset in -50:50 if offset != 0])
                sort!(phases)
            end
            results = [coupling_corrections(evolution_matrix(phase, modes; scheme, forcing, solver, filter=averaging),
                                           uncoupled_matrix(phase, modes; scheme, solver, filter=averaging), modes) for phase in phases]
            predictions = [frozen_predictions(phase, modes; scheme, solver, filter=averaging) for phase in phases]
            color = ("#666666", "#AA4499", "#0072B2", "#D55E00", "#009E73")[index] #src
            label = ("Exact, unfiltered", "Exact + F3", "RK3 + F3", "FB + F3", "FB + F1")[index] #src
            linestyle = (:solid, :dash, :solid, :dash, :dashdot)[index] #src
            for row in 1:3
                measured = [row == 1 ? result.direct_magnitude / result.damping : row == 2 ? result.indirect_magnitude / result.damping : result.growth / result.damping for result in results]
                predicted = [(row == 1 ? prediction.direct : row == 2 ? prediction.indirect : prediction.total) * 48 / (phase^4 * modes.χ[2]^2) - (row == 3) for (phase, prediction) in zip(phases, predictions)]
                lines!(axes[row], Float64.(phases), Float64.(measured); label, color=(color, 0.72), linestyle, linewidth=1.1) #src
                targets = [exp10.(range(log10(0.05), 0; length=9)) .* (1 + 0.035 * (index - 1)); #src
                           sort(unique([collect(1.4:0.65:20); collect(4.8:0.3:7.1)])) .+ 0.07 * (index - 1)] #src
                indices = unique([argmin(abs.(phases .- target)) for target in targets]) #src
                scatter!(axes[row], Float64.(phases[indices]), Float64.(predicted[indices]); color=:white, #src
                         strokecolor=color, strokewidth=0.7, markersize=3.4) #src
            end
            for (phase, result, prediction) in zip(phases, results, predictions)
                push!(rows, (; scheme, response=response_name, phase, direct=result.direct_magnitude, indirect=result.indirect_magnitude,
                              growth=result.growth, damping=result.damping, predicted_direct=prediction.direct, predicted_total=prediction.total))
            end
        end
        for row in 1:3 #src
            limits = ((-1.08, 0.22), (-0.55, 0.7), (-2.5, 0.5))[row] #src
            limits!(axes[row], 0, 20, limits...) #src
            hlines!(axes[row], [0]; color=:gray, linewidth=0.6) #src
            row == 3 && hlines!(axes[row], [-1]; color=:gray, linestyle=:dot, linewidth=0.7) #src
        end #src
        if column == 2 #src
            lines!(axes[3], [NaN], [NaN]; color=(:gray, 0.72), label=L"\mathrm{Matrix},\ \varepsilon=0.01") #src
            scatter!(axes[3], [NaN], [NaN]; color=:white, strokecolor=:black, strokewidth=0.7, markersize=4, label="Leading approximation") #src
            Legend(figure[4, 1:2], axes[3]; orientation=:horizontal, nbanks=3, labelsize=8, patchsize=(20, 8), rowgap=4, padding=(5, 5, 5, 5)) #src
        end #src
    end
    Label(figure[1:3, 0], L"\mathrm{Contribution}/(1-|\lambda_j|)"; rotation=π/2) #src
    save_results(output_directory, "damping_balance", figure, rows) #src
end

#md # ![](damping_balance.png)

# ## Forcing reconstruction and MRK4
#
# The stratification experiment uses μ₀ = 2 and exact, continuously forced barotropic integration. It compares frozen,
# final-stage quadratic, and progressive forcing with the leading analytical prediction. The MRK4 figure compares the
# direct correction and the full eigenvalue correction, scaled by θ₁⁴, using RK3 substeps and F3 at ε = 0.05.

setprecision(256) do
    reference_phase = big"2.0"
    figure = Figure(size=(518, 191), fontsize=8) #src
    axes = [Axis(figure[1, column], xscale=log10, yscale=log10, xlabel=L"\varepsilon", #src
                 ylabel=column == 1 ? "Absolute magnitude correction" : "", title=("Direct", "Barotropic indirect", "Combined correction")[column]) for column in 1:3] #src
    rows = NamedTuple[]
    insets = [Axis(figure[1, column], width=Relative(0.36), height=Relative(0.27), halign=0.95, valign=0.08, #src
                   xscale=log10, title="Frozen / prediction", titlesize=6.5, xticklabelsize=6.5, yticklabelsize=6.5, #src
                   xticks=[1e-3, 1e-1], limits=((1e-3, 1), (0.7, 1.05))) for column in 1:3] #src
    configurations = ((ExactSolver(), NoFilter(), ContinuousSampling()), (RungeKutta3(), F3Filter(), MidpointSampling()),
                      (ForwardBackward(), F3Filter(), MidpointSampling()))
    for (forcing_index, reconstruction) in enumerate((FrozenForcing, FinalStageQuadraticForcing, ProgressiveForcing)), (solver, averaging, sampling) in configurations
        forcing = reconstruction === FrozenForcing ? FrozenForcing() : reconstruction(; sampling)
        response_name = sampling isa ContinuousSampling ? "continuous" : solver isa RungeKutta3 ? "RK3_F3" : "FB_F3"
        filter = averaging_filter(averaging; float_type=BigFloat)
        stratifications = [big"0.64" / 4^index for index in 0:10]
        results = map(stratifications) do stratification
            modes = vertical_modes(stratification, 4)
            matrix = evolution_matrix(reference_phase, modes; forcing, solver, filter)
            reference = uncoupled_matrix(reference_phase, modes; solver, filter)
            result = coupling_corrections(matrix, reference, modes)
            prediction = frozen_predictions(reference_phase, modes; solver, filter)
            push!(rows, (; forcing, response=response_name, phase=reference_phase, stratification,
                          direct=result.direct_magnitude, indirect=result.indirect_magnitude, other=result.other_magnitude,
                          coupling=result.coupling, predicted_direct=prediction.direct, predicted_indirect=prediction.indirect))
            return result
        end
        if sampling isa ContinuousSampling
            color = ("#0072B2", "#D55E00", "#009E73")[forcing_index] #src
            label = ("Frozen", "Final-stage quadratic", "Progressive")[forcing_index] #src
            marker = (:circle, :rect, :utriangle)[forcing_index] #src
            for column in 1:3
                measured = [column == 1 ? result.direct_magnitude : column == 2 ? result.indirect_magnitude : result.coupling for result in results]
                scatterlines!(axes[column], Float64.(stratifications), Float64.(abs.(measured)); label, color, marker, #src
                              markercolor=:white, strokecolor=color, strokewidth=0.8, markersize=3, linewidth=0.9) #src
                if forcing in (FrozenForcing(), ProgressiveForcing()) #src
                    exponent = forcing == FrozenForcing() ? 2 : 3 #src
                    guide = reverse(stratifications)[3:5] #src
                    amplitude = abs(reverse(measured)[4]) * 8 / guide[2]^exponent #src
                    lines!(axes[column], Float64.(guide), Float64.(amplitude .* guide.^exponent); color=:gray, linestyle=:dot) #src
                    text!(axes[column], Float64(guide[end]), Float64(amplitude * guide[end]^exponent); #src
                          text=exponent == 2 ? L"\varepsilon^2" : L"\varepsilon^3", align=(:left, :bottom)) #src
                end #src
            end
            if forcing == FrozenForcing()
                predictions = [frozen_predictions(reference_phase, vertical_modes(stratification, 4); solver, filter) for stratification in stratifications]
                for column in 1:3
                    measured = [column == 1 ? result.direct_magnitude : column == 2 ? result.indirect_magnitude : result.coupling for result in results]
                    predicted = [column == 1 ? prediction.direct : column == 2 ? prediction.indirect : prediction.total for prediction in predictions]
                    ## Leading weak-stratification curves use χ₁² = 2ε²/π⁴.
                    for index in eachindex(stratifications)
                        modes = vertical_modes(stratifications[index], 4)
                        predicted[index] *= (2stratifications[index]^2 / BigFloat(π)^4) / modes.χ[2]^2
                    end
                    lines!(axes[column], Float64.(stratifications), Float64.(abs.(predicted)); color=:black, linestyle=:dash) #src
                    scatterlines!(insets[column], Float64.(stratifications), Float64.(measured ./ predicted); color="#0072B2", markercolor=:white, strokecolor="#0072B2", strokewidth=0.6, markersize=3) #src
                    hlines!(insets[column], [1]; color=:gray, linestyle=:dot) #src
                end
            end
        end
    end
    lines!(axes[1], [NaN], [NaN]; color=:black, linestyle=:dash, label="Analytical prediction (frozen)") #src
    Legend(figure[2, 1:3], axes[1]; orientation=:horizontal, nbanks=2, labelsize=7.5, patchsize=(20, 8), rowgap=2, padding=(4, 4, 4, 4)) #src
    save_results(output_directory, "forcing_stratification", figure, rows) #src
end

#md # ![](forcing_stratification.png)

# ### MRK4 coupling
#
# We normalize the corrections by θ₁⁴ and estimate the small-phase limit from two successively halved phases.

setprecision(256) do
    figure = Figure(size=(518, 230), fontsize=10) #src
    direct_axis = Axis(figure; bbox=BBox(62.2, 252.1, 69, 211.6), xlabel=L"\mu_0", #src
                       ylabel=L"\mathrm{Re}(\overline{\lambda_j}V_{jj}) / (|\lambda_j|\theta_1^4)", title="Direct", ylabelsize=11) #src
    total_axis = Axis(figure; bbox=BBox(320.4, 510.2, 69, 211.6), xlabel=L"\mu_0", #src
                      ylabel=L"(|\lambda_j^{\mathrm{full}}|-|\lambda_j|)/\theta_1^4", title="Combined correction", ylabelsize=11) #src
    for axis in (direct_axis, total_axis) #src
        axis.xscale = Makie.Symlog10(0.1; linscale=10/9) #src
        axis.xticks = ([0, 0.1, 1, 5, 10], ["0", "0.1", "1", "5", "10"]) #src
        axis.yticks = [-0.02, 0, 0.02] #src
        limits!(axis, 0, 10, -0.026, 0.034) #src
        hlines!(axis, [0]; color=:gray, linewidth=0.6) #src
        hlines!(axis, [-1 / 48, 1 / 48]; color=:gray, linestyle=:dot, linewidth=0.8) #src
        text!(axis, 0.045, 1 / 48; text=L"+1/48", color=:gray, fontsize=8.5, offset=(0, 3), align=(:center, :bottom)) #src
        text!(axis, 0.045, -1 / 48; text=L"-1/48", color=:gray, fontsize=8.5, offset=(0, 3), align=(:center, :bottom)) #src
    end #src
    hlines!(total_axis, [1 / 300]; color=:gray, linestyle=:dashdot, linewidth=0.8) #src
    text!(total_axis, 0.01, 1 / 300; text=L"1/300", color=:gray, fontsize=8.5, offset=(0, 3), align=(:left, :bottom)) #src
    modes = vertical_modes(big"0.05", 2)
    filter = averaging_filter(F3Filter(); float_type=BigFloat)
    rows = NamedTuple[]
    phases = samples == 9 ? exp10.(range(big"-2", big"1"; length=9)) :
             sort(unique([exp10.(range(big"-2", big"0"; length=18)); range(big"1.2", big"10"; length=45)]))
    for (index, forcing) in enumerate((FrozenForcing(), FinalStageQuadraticForcing(), ProgressiveForcing()))
        color = ("#0072B2", "#D55E00", "#009E73")[index] #src
        label = ("Frozen", "Final-stage quadratic", "Progressive")[index] #src
        results = map(phases) do phase
            matrix = evolution_matrix(phase, modes; scheme=ModifiedRungeKutta4(), forcing, solver=RungeKutta3(), filter)
            reference = uncoupled_matrix(phase, modes; scheme=ModifiedRungeKutta4(), solver=RungeKutta3(), filter)
            result = coupling_corrections(matrix, reference, modes)
            scale = (phase * modes.α[2])^4
            push!(rows, (; forcing, phase, direct=result.direct_magnitude / scale, indirect=result.indirect_magnitude / scale,
                          retained=(result.direct_magnitude + result.indirect_magnitude) / scale, net=result.coupling / scale, growth=result.growth / scale))
            return (result.direct_magnitude / scale, result.coupling / scale)
        end
        lines!(direct_axis, Float64.(phases), Float64.(first.(results)); label, color) #src
        lines!(total_axis, Float64.(phases), Float64.(last.(results)); label, color) #src
        limits = map((big"0.02", big"0.01")) do phase
            matrix = evolution_matrix(phase, modes; scheme=ModifiedRungeKutta4(), forcing, filter)
            reference = uncoupled_matrix(phase, modes; scheme=ModifiedRungeKutta4(), filter)
            result = coupling_corrections(matrix, reference, modes)
            return result.coupling / (phase * modes.α[2])^4
        end
        limit = (4limits[2] - limits[1]) / 3
        label = index == 1 ? L"\Gamma_j^{\mathrm{F}} \simeq %$(round(Float64(limit); digits=5))" : #src
                index == 2 ? L"\Gamma_j^{\mathrm{Q}} \simeq %$(round(Float64(limit); digits=5))" : #src
                             L"\Gamma_j^{\mathrm{P}} \simeq %$(round(Float64(limit) * 1e5; digits=2)) \times 10^{-5}" #src
        text!(total_axis, (0.40, 0.34, 0.30)[index], (0.14, 0.69, 0.32)[index]; text=label, color, space=:relative, fontsize=8.5) #src
    end
    predictions = [frozen_predictions(phase, modes; scheme=ModifiedRungeKutta4(), filter) for phase in phases]
    scales = (phases .* modes.α[2]).^4
    scatter!(direct_axis, Float64.(phases[1:3:end]), Float64.([prediction.direct for prediction in predictions][1:3:end] ./ scales[1:3:end]); #src
             color=:white, strokecolor="#0072B2", strokewidth=0.8, markersize=3.5, label="Frozen prediction") #src
    scatter!(total_axis, Float64.(phases[1:3:end]), Float64.([prediction.total for prediction in predictions][1:3:end] ./ scales[1:3:end]); #src
             color=:white, strokecolor="#0072B2", strokewidth=0.8, markersize=3.5) #src
    Legend(figure, direct_axis; bbox=BBox(60, 458, 3, 41), orientation=:horizontal, nbanks=2, labelsize=9, rowgap=2, colgap=15, patchsize=(20, 10), padding=(4, 4, 4, 4)) #src
    save_results(output_directory, "mrk4_forcing", figure, rows) #src
end

#md # ![](mrk4_forcing.png)

# ## Stability of the proposed configurations
#
# The full-spectrum scan compares all six configurations with eight and thirty-two baroclinic modes at ε = 0.05 and 0.15.
# Positive maximum growth indicates instability at that sampled phase. The final precision check repeats the weakly
# stratified forward–backward calculation at 128 and 256 bits to distinguish small growth from round-off.

setprecision(128) do
    figure = Figure(size=(540, 390), fontsize=10, figure_padding=(12, 16, 10, 14)) #src
    rows = NamedTuple[]
    names = ["AB2-SE", "WRK3-SE-F1", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM"]
    for (column, count) in enumerate((8, 32)), (row, stratification) in enumerate((big"0.05", big"0.15"))
        modes = vertical_modes(stratification, count)
        floating_modes = vertical_modes(Float64(stratification), count)
        axis = Axis(figure[row, column], xlabel=row == 2 ? L"\mu_0" : "", ylabel=column == 1 ? L"\max |\lambda| - 1" : "", #src
                    xticklabelsvisible=row == 2, yticklabelsvisible=column == 1, #src
                    xticks=0:10:40, titlegap=8, yticklabelspace=39.0, #src
                    xgridcolor=(:black, 0.08), ygridcolor=(:black, 0.08), #src
                    yscale=Makie.Symlog10(1e-13; linscale=0.35 * 10/9), #src
                    yticks=([-1e-2, -1e-6, -1e-10, 0, 1e-6], [L"-10^{-2}", L"-10^{-6}", L"-10^{-10}", L"0", L"10^{-6}"]), #src
                    title=L"\varepsilon=%$(Float64(stratification)),\ N_m=%$count") #src
        limits!(axis, 0, 44, -6e-2, 3e-4) #src
        phases = samples == 9 ? range(big"0.25", big"44"; length=9) : big"0.25":big"0.25":big"44"
        for (index, name) in enumerate(names)
            configuration = scheme_configuration(name)
            filter = averaging_filter(configuration.averaging; float_type=BigFloat)
            floating_filter = averaging_filter(configuration.averaging)
            growth = map(phases) do phase
                matrix = evolution_matrix(Float64(phase), floating_modes; configuration.scheme, configuration.forcing,
                                          configuration.solver, filter=floating_filter)
                value = maximum(abs, eigvals(matrix)) - 1
                bits = 53
                if abs(value) < 1e-12
                    matrix = evolution_matrix(phase, modes; configuration.scheme, configuration.forcing, configuration.solver, filter)
                    value = maximum(abs, eigvals(matrix)) - 1
                    bits = precision(BigFloat)
                end
                push!(rows, (; scheme=name, phase, stratification, modes=count, growth=value, precision=bits))
                return Float64(value)
            end
            lines!(axis, Float64.(phases), growth; label=name, #src
                   color=("#0000ff", "#ff0000", "#ff0000", "#008000", "#800080", "#ff0000")[index], #src
                   linestyle=(:solid, :dot, :solid, :dash, :solid, :dashdot)[index], linewidth=index == 2 ? 1.1 : 0.9) #src
            @info "Stability" scheme=name modes=count stratification=Float64(stratification) #src
        end
        vlines!(axis, Float64.([29 // 100 * sqrt(big"3"), sqrt(big"3"), 2sqrt(big"2.0")]) / Float64(modes.α[2]); color=:gray, linestyle=:dot) #src
        row == 2 && column == 2 && Legend(figure[3, 1:2], axis; nbanks=2, orientation=:horizontal, labelsize=9, patchsize=(25, 10), rowgap=5, colgap=18, padding=(8, 8, 6, 6)) #src
    end
    rowgap!(figure.layout, 1, 22) #src
    rowgap!(figure.layout, 2, 14) #src
    colgap!(figure.layout, 12) #src
    save_results(output_directory, "schemes", figure, rows) #src
end

let
    rows = NamedTuple[]
    for bits in (128, 256), count in (8, 32)
        setprecision(bits) do
            modes = vertical_modes(big"0.00002", count)
            filter = averaging_filter(F3Filter(); float_type=BigFloat)
            for scheme in (WickerSkamarock3(), ShuOsher3())
                phase = scheme == WickerSkamarock3() ? big"0.4" : big"0.35"
                forcing = scheme == ShuOsher3() ? SSPForcing() : FrozenForcing()
                matrix = evolution_matrix(phase, modes; scheme, forcing, solver=ForwardBackward(), filter)
                push!(rows, (; precision=bits, modes=count, scheme, phase, growth=maximum(abs, eigvals(matrix)) - 1))
            end
        end
    end
    open(joinpath(output_directory, "precision.csv"), "w") do stream #src
        println(stream, join(string.(keys(first(rows))), ',')) #src
        for row in rows #src
            println(stream, join(values(row), ',')) #src
        end #src
    end #src
end

#md # ![](schemes.png)

# ## Global temporal accuracy
#
# Each configuration advances the exact initial state to k√(gH)T = 64. AB2 receives an exact previous state. The error
# is the operator 2-norm in the energy-scaled coordinates [u; h/α]. The timestep is adjusted to reach the final time exactly.

setprecision(256) do
    figure = Figure(size=(461, 220), fontsize=9) #src
    axis = Axis(figure; bbox=BBox(50.7, 336.5, 39.6, 213.4), xlabel=L"\mu_0", ylabel="Global error", xscale=log10, yscale=log10, xlabelsize=10, ylabelsize=10) #src
    limits!(axis, 8e-4, 0.6, 1e-10, 4) #src
    axis.yticks = (10.0 .^ (-10:2:0), [L"10^{%$exponent}" for exponent in -10:2:0]) #src
    axis.xgridvisible = false #src
    rows = NamedTuple[]
    modes = vertical_modes(big"0.05", 2)
    phases = exp10.(range(big"-3", log10(big"0.5"); length=min(samples, 25)))
    for (index, name) in enumerate(["AB2-SE", "WRK3-SE-F1", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM"])
        configuration = scheme_configuration(name)
        filter = averaging_filter(configuration.averaging; float_type=BigFloat)
        results = [global_error(phase, modes; configuration.scheme, configuration.forcing, configuration.solver, filter) for phase in phases]
        for result in results
            push!(rows, (; scheme=name, phase=result.phase, error=result.error))
        end
        color = ("#0000ff", "#ff0000", "#ff0000", "#008000", "#800080", "#ff0000")[index] #src
        lines!(axis, Float64.([result.phase for result in results]), [result.error for result in results]; label=name, color, #src
               linestyle=(:solid, :dot, :solid, :dash, :solid, :dashdot)[index]) #src
        if index in (2, 5) #src
            scatter!(axis, Float64.([result.phase for result in results][2:4:end]), [result.error for result in results][2:4:end]; #src
                     marker=index == 2 ? :rect : :diamond, color=:white, strokecolor=color, strokewidth=0.9, markersize=4, label=name) #src
        end #src
    end
    for exponent in 1:2 #src
        guide = [exponent == 1 ? 0.0015 : 0.004, 0.025] #src
        reference = (exponent == 1 ? 0.18 : 1e-5) .* guide.^exponent #src
        lines!(axis, guide, reference; color=:gray, linestyle=:dash, linewidth=0.8) #src
        text!(axis, guide[end], reference[end]; offset=(3, 4), align=(:left, :bottom), fontsize=8, color=:gray, #src
              text=exponent == 1 ? L"\mathrm{First\ order}:\ \mu_0" : L"\mathrm{Second\ order}:\ \mu_0^2") #src
    end #src
    Legend(figure, axis; bbox=BBox(348, 461, 45, 205), merge=true, labelsize=9, rowgap=6, patchsize=(20, 9), padding=(7, 7, 7, 7)) #src
    save_results(output_directory, "order_schemes", figure, rows) #src
end

#md # ![](order_schemes.png)

# Each figure is saved as PDF and PNG, with its numerical values in a companion CSV. The package tests independently
# check the direct stage updates, modal identities, filter moments, analytical transfers, and small-step convergence.
