# # Time discretization of free-surface ocean models
#
# This example reproduces the amplification, filtering, and coupled modal analyses of Sections 2–4. The numerical implementation is provided by TimestepperTheory. Each experiment shows its computations and passes the results to a plotting function in `time_discretization_figures.jl`.

using TimestepperTheory
using LinearAlgebra

include("time_discretization_figures.jl") #src

output_directory = joinpath(@__DIR__, "..", "output")
samples = "--quick" in ARGS ? 9 : 61
setprecision(256)

# ## Amplification factors of the outer schemes
#
# For a Fourier mode of unit-speed advection, the exact timestep multiplies the state by exp(−iθ), where θ = kΔt.
# The magnitude departure |λ| − 1 measures numerical damping or growth. WRK3 and SRK3 share the same cubic polynomial;
# RK2 and the fourth-order Taylor polynomial provide references. The AB2 curves use the physical root, with ξ = 0 and 0.1.

# ```math
# \theta = k\Delta t, \qquad \lambda = R(-i\theta), \qquad
# R_{\mathrm{MRK4}}(z) = 1 + z + \frac{z^2}{2} + \frac{z^3}{6} + \frac{a z^4}{6}.
# ```

let
    rows   = NamedTuple[]
    panels = NamedTuple[]
    phases = BigFloat.(0:300) ./ 100
    for ξ in (0, 1 // 10)
        growth = [abs(amplification_factor(θ, AdamsBashforth2(); ξ)) - 1 for θ in phases]
        push!(panels, (; phases, growth))
        append!(rows, [(; scheme="AB2", parameter=Float64(ξ), phase=θ, growth=value) for (θ, value) in zip(phases, growth)])
    end
    for (index, scheme) in enumerate((RungeKutta2(), WickerSkamarock3(), ModifiedRungeKutta4(a=1 // 4), ModifiedRungeKutta4()))
        growth = [abs(amplification_factor(θ, scheme)) - 1 for θ in phases]
        name   = ("RK2", "RK3", "RK4", "MRK4")[index]
        push!(panels, (; phases, growth))
        append!(rows, [(; scheme=name, parameter=scheme isa ModifiedRungeKutta4 ? Float64(scheme.a) : NaN, phase=θ, growth=value) for (θ, value) in zip(phases, growth)])
    end
    plot_amplification_factor(panels, rows, output_directory)
end

#md # ![](amplification_factor.png)

# ## Barotropic substeps and averaging filters
#
# A barotropic solve averages the substep evolution with weights w. The three panels compare exact, forward–backward,
# and RK3 substeps, using the F1, F2, and F3 filters and 64 nominal substeps. The grey curve is the implicit solver.

let
    phases          = sort(unique([0; exp10.(range(-2, 0; length=180)); range(1, 16; length=901)]))
    implicit_growth = @. inv(sqrt(1 + phases^2)) - 1
    rows            = NamedTuple[]
    panels          = NamedTuple[]
    for (column, solver) in enumerate((ExactSolver(), ForwardBackward(), RungeKutta3()))
        for (index, averaging) in enumerate((F1Filter(), F2Filter(), F3Filter()))
            filter = averaging_filter(averaging, 64)
            growth = map(phases) do μ₀
                response = barotropic_response(μ₀; solver, filter)
                return sum(abs, eigvals(response.W)) / 2 - 1
            end
            push!(panels, (; phases, growth))
            name = ("Exact", "FB", "RK3")[column]
            append!(rows, [(; solver=name, filter="F$index", phase=μ₀, growth=value) for (μ₀, value) in zip(phases, growth)])
        end
    end
    push!(panels, (; phases, growth=implicit_growth))
    append!(rows, [(; solver="Implicit", filter="none", phase=μ₀, growth=value) for (μ₀, value) in zip(phases, implicit_growth)])
    plot_barotropic_filter_damping(panels, rows, output_directory)
end

#md # ![](barotropic-filter-damping.png)

# The transport weights follow by summing the discrete continuity updates: each weight is the fractional spacing
# multiplied by the sum of all subsequent output weights, including the current one. Both sets of weights integrate
# to one. The centered moments M₂, M₃, and M₄ describe the filter's low-frequency damping and phase error.

# ```math
# w_s^\dagger = \Delta\tau\sum_{r=s}^{S}w_r, \qquad
# M_p = \sum_s w_s(\tau_s-1)^p.
# ```

let
    rows   = NamedTuple[]
    panels = NamedTuple[]
    for (column, averaging) in enumerate((F1Filter(), F2Filter(), F3Filter()))
        filter = averaging_filter(averaging, 64)
        w, Δτ  = filter.w, filter.Δτ
        τ      = eachindex(w) .* Δτ
        w⁺     = Δτ .* reverse(cumsum(reverse(w)))
        Mₚ     = [sum(w .* (τ .- 1).^p) for p in 2:4]

        push!(panels, (; w, Δτ, τ, w⁺, Mₚ))
        append!(rows, [(; filter="F$column", time=τ[s], spacing=Δτ, weight=w[s], transport_weight=w⁺[s])
                       for s in eachindex(w)])
    end
    plot_averaging_weights(panels, rows, output_directory)
end

#md # ![](averaging-weights.png)

# ## Vertical modes and the modal state
#
# The dimensionless stratification is ε = N²H/g. Each root of x tan(x) = ε determines a normalized wave speed α,
# a surface value ψ, and a depth average χ. The requested count excludes the barotropic mode, which is always retained.

# ```math
# x_q\tan x_q = \varepsilon, \qquad \alpha_q = \frac{\sqrt{\varepsilon}}{x_q}, \qquad
# \psi_q = A_q\cos x_q, \qquad \chi_q = \psi_q\alpha_q^2.
# ```

ε     = big"0.05"
modes = vertical_modes(ε, 2)

modes.α

# The state is ordered as y = [u₀, …, uQ, h₀, …, hQ], with velocity divided by √(gH) and pressure divided by gH.
# The depth average is obtained by multiplying the corresponding amplitude vector by χᵀ. At finite mode count,
# χᵀχ is slightly less than one. The modal weights retain this truncation error; the SRK3 blended matrices use χᵀχ ≈ 1.

omitted_weight = 1 - sum(abs2, modes.χ)

# ## Barotropic integration and polynomial forcing
#
# The barotropic Courant number μ₀ = k√(gH)Δt uses the reference speed √(gH). The averaging filter
# determines the returned velocity and continuity transport. F3 has vanishing second and third centered moments.

μ₀       = big"2.0"
filter   = averaging_filter(F3Filter(); float_type=BigFloat)
response = barotropic_response(μ₀; solver=RungeKutta3(), filter)

response.W

# The forcing outputs correspond to 1, s, and s², where s is time divided by the stage duration. Multiplying these
# outputs by the reconstructed polynomial coefficients gives the forced velocity and continuity transport.

response.ūₚ, response.ūₚ⁺

# The continuity transport follows from the same discrete surface update, including the averaging filter.

continuity_residual = norm(response.W[2, :] - [0, 1] + im * μ₀ * response.W̃)

# ## Stage matrices and the complete timestep
#
# For WRK3-SE, the stages use progressive reconstruction of the slow forcing. The returned arrays contain
# A[m] = Aₘ, B[m] = Bₘ, and M[m] = Mₘ. The initial matrix M₀ = I is used internally. Contributions from earlier stages enter through Aₘ.

# ```math
# \boldsymbol y^m = \mathbf A_m\boldsymbol y^n + \mathbf B_m\boldsymbol y^{m-1}, \qquad
# \mathbf M_0 = \mathbf I, \qquad \mathbf M_m = \mathbf A_m + \mathbf B_m\mathbf M_{m-1}.
# ```

stages = stage_matrices(μ₀, modes; scheme=WickerSkamarock3(), forcing=ProgressiveForcing(), solver=RungeKutta3(), filter)
Aₘ = stages.A[end]
Bₘ = stages.B[end]
Mₘ = stages.M[end]

Mₘ

# SRK3 uses its blended stage update with frozen forcing. MRK4 uses progressive reconstruction, including a
# least-squares quadratic at the fourth stage. The implicit WRK3 scheme uses an implicit surface solve. AB2 uses an
# augmented state [yⁿ; yⁿ⁻¹], so its complete matrix has twice the dimension of a Runge–Kutta timestep matrix.

configurations = ("AB2-SE", "WRK3-SE-F1", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM")

for name in configurations
    configuration = scheme_configuration(name)
    averaging     = averaging_filter(configuration.averaging; float_type=BigFloat)
    M             = evolution_matrix(μ₀, modes; configuration.scheme, configuration.forcing, configuration.solver, filter=averaging)
    println(name, ": matrix size ", size(M), ", maximum growth ", Float64(maximum(abs, eigvals(M)) - 1))
end

# ## Direct and indirect coupling corrections
#
# The uncoupled matrix advances each mode independently. Projecting the difference between the complete and uncoupled
# matrices onto their normalized left and right eigenvectors gives the coupling coefficients Vⱼₗ. These determine the
# direct correction within the selected mode and propagation direction, and the indirect correction through other modes.

# ```math
# V_{jl} = \ell_j(\mathbf M-\mathbf M_u)r_l, \qquad
# \delta\lambda_j^{(1)} = V_{jj}, \qquad
# \delta\lambda_j^{(2)} = \sum_{l\ne j}\frac{V_{jl}V_{lj}}{\lambda_j-\lambda_l}.
# ```
#
# The magnitude correction is obtained by projecting the eigenvalue shift onto the radial direction:
#
# ```math
# |\lambda_j+\delta\lambda_j|-|\lambda_j|
# \simeq \frac{\operatorname{Re}(\overline{\lambda_j}\,\delta\lambda_j)}{|\lambda_j|}.
# ```

Mᵤ = uncoupled_matrix(μ₀, modes; scheme=WickerSkamarock3(), solver=RungeKutta3(), filter)
corrections = coupling_corrections(Mₘ, Mᵤ, modes)

corrections.direct_magnitude, corrections.indirect_magnitude, corrections.other_magnitude

# The first two figures compare the perturbation expansion with the full eigenvalue and separate the damping balance
# for exact, RK3, and forward–backward barotropic integration. These calculations use 128-bit arithmetic and eight
# baroclinic modes. The decomposition includes both stratifications shown in the manuscript figure.

setprecision(128) do
    filter = averaging_filter(F3Filter(); float_type=BigFloat)
    phases = samples == 9 ? BigFloat.(range(0.1, 16; length=9)) : BigFloat.(1:160) ./ 10
    rows   = NamedTuple[]
    panels = NamedTuple[]
    for ε in (big"0.0005", big"0.05"), scheme in (WickerSkamarock3(), ShuOsher3())
        forcing     = FrozenForcing()
        modes       = vertical_modes(ε, 8)
        corrections = [coupling_corrections(evolution_matrix(μ₀, modes; scheme, forcing, filter),
                                            uncoupled_matrix(μ₀, modes; scheme, filter), modes) for μ₀ in phases]
        predictions = [frozen_predictions(μ₀, modes; scheme, filter) for μ₀ in phases]
        complete    = [-1 + (result.direct_magnitude + result.indirect_magnitude + result.other_magnitude) / result.damping
                     for result in corrections]
        retained  = [-1 + (result.direct_magnitude + result.indirect_magnitude) / result.damping for result in corrections]
        predicted = [-1 + real(result.reference_eigenvalue) / abs(result.reference_eigenvalue) * prediction.total / result.damping
                     for (prediction, result) in zip(predictions, corrections)]
        push!(panels, (; phases, corrections, complete, retained, predicted))
        for (index, μ₀) in enumerate(phases)
            result = corrections[index]
            push!(rows, (; scheme, phase=μ₀, stratification=ε, direct=result.direct_magnitude, indirect=result.indirect_magnitude,
                          other=result.other_magnitude, growth=result.growth, damping=result.damping, prediction=predicted[index], overlap=result.overlap))
        end
    end
    plot_coupling_decomposition(panels, rows, output_directory)
end

#md # ![](coupling_decomposition.png)

# ### Damping balance
#
# We repeat the decomposition with each barotropic solver and filter, and compare its direct and indirect contributions.

setprecision(128) do
    rows   = NamedTuple[]
    panels = NamedTuple[]
    for scheme in (WickerSkamarock3(), ShuOsher3())
        forcing        = FrozenForcing()
        configurations = ((ExactSolver(), NoFilter()), (ExactSolver(), F3Filter()), (RungeKutta3(), F3Filter()),
                          (ForwardBackward(), F3Filter()), (ForwardBackward(), F1Filter()))
        for (index, (solver, averaging)) in enumerate(configurations)
            response_name = ("exact", "exact_F3", "RK3_F3", "FB_F3", "FB_F1")[index]
            averaging     = averaging_filter(averaging; float_type=BigFloat)
            modes         = vertical_modes(big"0.01", 8)
            phases        = samples == 9 ? BigFloat.(range(0.05, 20; length=9)) : BigFloat.(sort(unique([collect(0.05:0.025:1); collect(1:0.1:20)])))
            if response_name == "exact" && samples != 9
                append!(phases, [2BigFloat(π)*resonance + BigFloat(offset)/500 for resonance in 1:3 for offset in -50:50 if offset != 0])
                sort!(phases)
            end
            results = [coupling_corrections(evolution_matrix(μ₀, modes; scheme, forcing, solver, filter=averaging),
                                           uncoupled_matrix(μ₀, modes; scheme, solver, filter=averaging), modes) for μ₀ in phases]
            predictions = [frozen_predictions(μ₀, modes; scheme, solver, filter=averaging) for μ₀ in phases]
            measured    = [[result.direct_magnitude / result.damping for result in results],
                        [result.indirect_magnitude / result.damping for result in results],
                        [result.growth / result.damping for result in results]]
            predicted = [[prediction.direct * 48 / (μ₀^4 * modes.χ[2]^2) for (μ₀, prediction) in zip(phases, predictions)],
                         [prediction.indirect * 48 / (μ₀^4 * modes.χ[2]^2) for (μ₀, prediction) in zip(phases, predictions)],
                         [prediction.total * 48 / (μ₀^4 * modes.χ[2]^2) - 1 for (μ₀, prediction) in zip(phases, predictions)]]
            push!(panels, (; phases, measured, predicted))
            for (μ₀, result, prediction) in zip(phases, results, predictions)
                push!(rows, (; scheme, response=response_name, phase=μ₀, direct=result.direct_magnitude, indirect=result.indirect_magnitude,
                              growth=result.growth, damping=result.damping, predicted_direct=prediction.direct, predicted_total=prediction.total))
            end
        end
    end
    plot_damping_balance(panels, rows, output_directory)
end

#md # ![](damping_balance.png)

# ## Forcing reconstruction and MRK4
#
# The slow forcing is reconstructed from the previously computed stage pressures. Linearity gives the forced outputs:
#
# ```math
# F(t') = \sum_p F_p t'^p, \qquad
# \bar u_{\mathrm{forced}}^m = \sum_p F_p\bar u_p, \qquad
# \bar u_{\mathrm{forced}}^\dagger = \sum_p F_p\bar u_p^\dagger.
# ```
#
# The stratification experiment uses μ₀ = 2 and exact, continuously forced barotropic integration. It compares frozen
# and progressive forcing. Open circles show the frozen-forcing analytical prediction; the lower panels compare the
# computed frozen correction with this prediction. The dotted guides show ε² and ε³ scaling.

setprecision(256) do
    μ₀             = big"2.0"
    rows           = NamedTuple[]
    panels         = NamedTuple[]
    configurations = ((ExactSolver(), NoFilter(), ContinuousSampling()), (RungeKutta3(), F3Filter(), MidpointSampling()),
                      (ForwardBackward(), F3Filter(), MidpointSampling()))
    for reconstruction in (FrozenForcing, ProgressiveForcing), (solver, averaging, sampling) in configurations
        forcing         = reconstruction === FrozenForcing ? FrozenForcing() : reconstruction(; sampling)
        response_name   = sampling isa ContinuousSampling ? "continuous" : solver isa RungeKutta3 ? "RK3_F3" : "FB_F3"
        filter          = averaging_filter(averaging; float_type=BigFloat)
        stratifications = [big"0.64" / 4^index for index in 0:10]
        results         = map(stratifications) do ε
            modes      = vertical_modes(ε, 4)
            M          = evolution_matrix(μ₀, modes; forcing, solver, filter)
            Mᵤ         = uncoupled_matrix(μ₀, modes; solver, filter)
            result     = coupling_corrections(M, Mᵤ, modes)
            prediction = forcing isa FrozenForcing ? frozen_predictions(μ₀, modes; solver, filter) :
                         (; direct=BigFloat(NaN), indirect=BigFloat(NaN))
            push!(rows, (; forcing, response=response_name, phase=μ₀, stratification=ε,
                          direct=result.direct_magnitude, indirect=result.indirect_magnitude, other=result.other_magnitude,
                          coupling=result.coupling, predicted_direct=prediction.direct, predicted_indirect=prediction.indirect))
            return result
        end
        if sampling isa ContinuousSampling
            measured  = [[result.direct_magnitude for result in results],
                        [result.indirect_magnitude for result in results], [result.coupling for result in results]]
            predicted = Vector{BigFloat}[]
            if forcing isa FrozenForcing
                predictions = [frozen_predictions(μ₀, vertical_modes(ε, 4); solver, filter) for ε in stratifications]
                predicted   = [[prediction.direct for prediction in predictions],
                             [prediction.indirect for prediction in predictions], [prediction.total for prediction in predictions]]
                for values in predicted, index in eachindex(stratifications)
                    modes = vertical_modes(stratifications[index], 4)
                    values[index] *= (2stratifications[index]^2 / BigFloat(π)^4) / modes.χ[2]^2
                end
            end
            push!(panels, (; forcing, stratifications, measured, predicted))
        end
    end
    plot_forcing_stratification(panels, rows, output_directory)
end

#md # ![](forcing_stratification.png)

# ### MRK4 coupling
#
# We compare frozen and progressive forcing at a=1/4, retaining 24 baroclinic modes. The insets resolve the small-step
# residual with 192-bit arithmetic. Both this example and the standalone figure use the same experiment.

include("mrk4_forcing_prediction.jl")
mrk4_forcing_experiment(output_directory; samples)

#md # ![](mrk4_forcing.png)

# ## Stability of the proposed configurations
#
# The full-spectrum scan compares all six configurations with eight and thirty-two baroclinic modes at ε = 0.05 and 0.15.
# Positive maximum growth indicates instability at that sampled μ₀. The final precision check repeats the weakly
# stratified forward–backward calculation at 128 and 256 bits to distinguish small growth from round-off.

setprecision(128) do
    rows   = NamedTuple[]
    panels = NamedTuple[]
    names  = ["AB2-SE", "WRK3-SE-F1", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM"]
    for count in (8, 32), ε in (big"0.05", big"0.15")
        modes          = vertical_modes(ε, count)
        floating_modes = vertical_modes(Float64(ε), count)
        phases         = samples == 9 ? range(big"0.25", big"44"; length=9) : big"0.25":big"0.25":big"44"
        curves         = Vector{Float64}[]
        for name in names
            configuration   = scheme_configuration(name)
            filter          = averaging_filter(configuration.averaging; float_type=BigFloat)
            floating_filter = averaging_filter(configuration.averaging)
            growth          = map(phases) do μ₀
                M = evolution_matrix(Float64(μ₀), floating_modes; configuration.scheme, configuration.forcing,
                                          configuration.solver, filter=floating_filter)
                value = maximum(abs, eigvals(M)) - 1
                bits  = 53
                if abs(value) < 1e-12
                    M     = evolution_matrix(μ₀, modes; configuration.scheme, configuration.forcing, configuration.solver, filter)
                    value = maximum(abs, eigvals(M)) - 1
                    bits  = precision(BigFloat)
                end
                push!(rows, (; scheme=name, phase=μ₀, stratification=ε, modes=count, growth=value, precision=bits))
                return Float64(value)
            end
            push!(curves, growth)
        end
        thresholds = Float64.([29 // 100 * sqrt(big"3"), sqrt(big"3"), 2sqrt(big"2.0")]) / Float64(modes.α[2])
        push!(panels, (; phases, growth=curves, thresholds))
    end
    plot_schemes(panels, rows, output_directory)
end

let
    rows = NamedTuple[]
    for bits in (128, 256), count in (8, 32)
        setprecision(bits) do
            modes  = vertical_modes(big"0.00002", count)
            filter = averaging_filter(F3Filter(); float_type=BigFloat)
            for scheme in (WickerSkamarock3(), ShuOsher3())
                μ₀      = scheme == WickerSkamarock3() ? big"0.4" : big"0.35"
                forcing = FrozenForcing()
                M       = evolution_matrix(μ₀, modes; scheme, forcing, solver=ForwardBackward(), filter)
                push!(rows, (; precision=bits, modes=count, scheme, phase=μ₀, growth=maximum(abs, eigvals(M)) - 1))
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
    rows   = NamedTuple[]
    panels = Vector{NamedTuple}[]
    modes  = vertical_modes(big"0.05", 2)
    phases = exp10.(range(big"-3", log10(big"0.5"); length=min(samples, 25)))
    for name in ("AB2-SE", "WRK3-SE-F1", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM")
        configuration = scheme_configuration(name)
        filter        = averaging_filter(configuration.averaging; float_type=BigFloat)
        results       = [global_error(μ₀, modes; configuration.scheme, configuration.forcing, configuration.solver, filter) for μ₀ in phases]
        push!(panels, results)
        for result in results
            push!(rows, (; scheme=name, phase=result.phase, error=result.error))
        end
    end
    plot_order_schemes(panels, rows, output_directory)
end

#md # ![](order_schemes.png)

# Each figure is saved as PDF and PNG, with its numerical values in a companion CSV. The package tests independently
# check the direct stage updates, modal identities, filter moments, analytical transfers, and small-step convergence.
