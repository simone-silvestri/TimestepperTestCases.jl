using LinearAlgebra
using TimestepperTheory

function mrk4_forcing_experiment(output_directory; samples=61)
    BLAS.set_num_threads(1)
    return setprecision(192) do
        scheme   = ModifiedRungeKutta4(a=1 // 4)
        sampling = MidpointSampling()
        panels   = NamedTuple[]
        rows     = NamedTuple[]
        phases = samples == 9 ? BigFloat.(range(1, 30; length=9)) : collect(big"0.1":big"0.1":big"30")
        phases = sort(unique([exp10.(range(big"-3", big"0"; length=25)); phases]))
        prediction_phases = BigFloat.(["0.001", "0.01", "0.05", "0.15", "0.3", "0.5", "0.8", "1.5", "2.4", "3.2",
                                     "4", "4.6", "5", "5.4", "5.8", "6.2", "6.7", "7.4", "8.3", "9.5", "11", "13", "16", "20", "25", "29.5"])
        for ε in (big"0.0005", big"0.05")
            modes  = vertical_modes(ε, 24)
            filter = averaging_filter(F3Filter(); float_type=BigFloat)
            floating_modes  = vertical_modes(Float64(ε), 24)
            floating_filter = averaging_filter(F3Filter())
            results = map(phases) do μ₀
                # Resolve the small-step residual before dividing by θ₁⁴.
                high_precision = μ₀ <= 1
                current_modes  = high_precision ? modes : floating_modes
                current_filter = high_precision ? filter : floating_filter
                phase = high_precision ? μ₀ : Float64(μ₀)
                Mᵤ = uncoupled_matrix(phase, current_modes; scheme, solver=RungeKutta3(), filter=current_filter)
                corrections = map((FrozenForcing(), ProgressiveForcing(; sampling))) do forcing
                    M = evolution_matrix(phase, current_modes; scheme, forcing, solver=RungeKutta3(), filter=current_filter)
                    coupling_corrections(M, Mᵤ, current_modes)
                end
                scale = (μ₀ * modes.α[2])^4
                values = Float64.((corrections[1].coupling, corrections[2].coupling, corrections[1].damping) ./ scale)
                push!(rows, (; stratification=Float64(ε), phase=Float64(μ₀), kind="computed",
                              frozen=values[1], progressive=values[2], damping=values[3]))
                return values
            end
            predictions = map(prediction_phases) do μ₀
                frozen = frozen_predictions(μ₀, modes; scheme, solver=RungeKutta3(), filter)
                scale = (μ₀ * modes.α[2])^4
                value = Float64(frozen.total / scale)
                push!(rows, (; stratification=Float64(ε), phase=Float64(μ₀), kind="prediction",
                              frozen=value, progressive=NaN, damping=NaN))
                return value
            end
            θ = phases .* modes.α[2]
            damping = [(1 .- abs.(amplification_factor.(θ, Ref(ModifiedRungeKutta4(; a))))) ./ θ.^4 for a in (0, 23 // 100)]
            push!(panels, (; ε, phases, results, prediction_phases, predictions, damping))
            @info "Computed MRK4 coupling" ε=Float64(ε)
        end
        plot_mrk4_forcing(panels, rows, output_directory)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    include("time_discretization_figures.jl")
    output_directory = isempty(ARGS) ? joinpath(@__DIR__, "..", "output") : first(ARGS)
    mrk4_forcing_experiment(output_directory)
end
