function direct_step(phase, modes, initial; scheme=WickerSkamarock3(), forcing=FinalStageQuadraticForcing(), solver=RungeKutta3(),
                     filter=averaging_filter(F3Filter(); float_type=typeof(phase)), sampling=MidpointSampling(), a=23 // 100)
    count = length(modes.α)
    velocity = initial[1:count]
    pressure = initial[count+1:2count]
    fractions = scheme == ModifiedRungeKutta4() ? [a, 1 // 3, 1 // 2, 1] : scheme == ShuOsher3() ? [1, 1, 1] : [1 // 3, 1 // 2, 1]
    pressures = [copy(pressure)]
    initial_velocity = sum(modes.χ .* velocity)
    initial_surface = sum(modes.ψ .* pressure)
    for stage in eachindex(fractions)
        stage_phase = phase * fractions[stage]
        if solver == ImplicitSolver()
            surface = sum(modes.ψ .* pressure)
            predictor = initial[1:count] - im * stage_phase * (pressure - surface * modes.χ)
            surface = (initial_surface - im * stage_phase * sum(modes.χ .* predictor)) / (1 + stage_phase^2)
            next_velocity = predictor - im * stage_phase * surface * modes.χ
            transport = sum(modes.χ .* next_velocity)
            corrected_velocity = velocity + modes.χ * (transport - sum(modes.χ .* velocity))
            next_pressure = initial[count+1:2count] - im * stage_phase * modes.α.^2 .* corrected_velocity
        else
            response = barotropic_response(stage_phase; solver, filter, sampling)
            samples = [sum(modes.Π .* sample) for sample in pressures]
            coefficients = zeros(eltype(initial), 3)
            if forcing == SSPForcing() && stage == 3
                coefficients[1] = samples[1] / 6 + samples[2] / 6 + 2samples[3] / 3
            elseif forcing in (FrozenForcing(), SSPForcing()) || stage == 1 || (forcing == FinalStageQuadraticForcing() && stage < length(fractions))
                coefficients[1] = samples[end]
            else
                nodes = [0; fractions[1:stage-1]] ./ fractions[stage]
                degree = min(stage - 1, 2)
                vandermonde = [typeof(phase)(node)^power for node in nodes, power in 0:degree]
                coefficients[1:degree+1] = vandermonde \ samples
            end
            barotropic_velocity = response.W[1, 1] * initial_velocity + response.W[1, 2] * initial_surface + sum(response.ūₚ .* coefficients)
            transport = response.W̃[1] * initial_velocity + response.W̃[2] * initial_surface + sum(response.ūₚ⁺ .* coefficients)
            corrected_velocity = velocity + modes.χ * (transport - sum(modes.χ .* velocity))
            if scheme == ShuOsher3()
                blending = (1, 1 // 4, 2 // 3)[stage]
                increment = -im * stage_phase * pressure
                increment += modes.χ * (barotropic_velocity - initial_velocity - sum(modes.χ .* increment))
                next_velocity = (1 - blending) * initial[1:count] + blending * (velocity + increment)
                next_pressure = (1 - blending) * initial[count+1:2count] + blending * (pressure - im * stage_phase * modes.α.^2 .* corrected_velocity)
            else
                next_velocity = initial[1:count] - im * stage_phase * pressure
                next_pressure = initial[count+1:2count] - im * stage_phase * modes.α.^2 .* corrected_velocity
                next_velocity += modes.χ * (barotropic_velocity - sum(modes.χ .* next_velocity))
            end
        end
        velocity, pressure = next_velocity, next_pressure
        push!(pressures, pressure)
    end
    return [velocity; pressure]
end

function direct_adams_bashforth(phase, modes, history, filter)
    count = length(modes.α)
    initial = history[1:2count]
    previous = history[2count+1:end]
    extrapolated = (8initial - 3previous) / 5
    response = barotropic_response(phase; solver=ForwardBackward(), filter)
    initial_velocity = sum(modes.χ .* initial[1:count])
    initial_surface = sum(modes.ψ .* initial[count+1:end])
    forcing = sum(modes.Π .* extrapolated[count+1:end])
    barotropic_velocity = response.W[1, 1] * initial_velocity + response.W[1, 2] * initial_surface + response.ūₚ[1] * forcing
    transport = response.W̃[1] * initial_velocity + response.W̃[2] * initial_surface + response.ūₚ⁺[1] * forcing
    velocity = initial[1:count] - im * phase * extrapolated[count+1:end]
    corrected_velocity = extrapolated[1:count] + modes.χ * (transport - sum(modes.χ .* extrapolated[1:count]))
    pressure = initial[count+1:end] - im * phase * modes.α.^2 .* corrected_velocity
    velocity += modes.χ * (barotropic_velocity - sum(modes.χ .* velocity))
    return [velocity; pressure; initial]
end
