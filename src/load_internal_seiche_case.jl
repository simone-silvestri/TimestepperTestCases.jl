"""
    load_internal_seiche(folder, timestepper, free_surface, μ₀)

Load an internal seiche run and reduce it to the two numbers the case exists to measure.

$(SIGNATURES)

# Arguments
- `folder`: Directory containing the simulation output files
- `timestepper`: Timestepper name string (e.g., `"SplitRungeKutta3"`)
- `free_surface`: Configuration name string (e.g., `"rk3_frozen"`)
- `μ₀`: Nominal barotropic Courant number the run was performed at

# Returns
- Dictionary containing:
  - `:KE`, `:APE`, `:E`: kinetic, available potential and total perturbation energy time series
  - `:iterations`, `:times`: the sampling iterations and times
  - `:λ`: measured damping rate of the mode amplitude, per time step
  - `:λ_err`: standard error of `λ` from the least-squares fit
  - `:ω`, `:μ₀_eff`: measured oscillation frequency and the *effective* barotropic Courant number
  - `:μ₀`: the nominal barotropic Courant number requested

The energy of a freely-decaying normal mode falls as `E(n) = E(0) exp(2λn)`, so `λ` is read off the
slope of `log E` against iteration and halved. Because the configuration carries no physical
dissipation, `λ` is the amplification factor of the coupled barotropic--baroclinic scheme.

`μ₀_eff` matters. A second-order centered discretization propagates the mode at the *discrete* phase
speed, whose effective wavenumber is `sin(kΔx)/Δx` rather than `k`, so the Courant number the scheme
actually experiences is smaller than the nominal one by up to a factor of `2/π` at the grid scale. The
effective value is recovered from the measured oscillation frequency of the mode, and is the abscissa
against which the analytic predictions should be plotted.
"""
function load_internal_seiche(folder, timestepper, free_surface, μ₀)
    tag  = replace(string(μ₀), "." => "p")
    path = folder * "internal_seiche_" * timestepper * "_" * free_surface * "_mu" * tag * "_energy.jld2"
    case = Dict()

    KE  = FieldTimeSeries(path, "KE")
    APE = FieldTimeSeries(path, "APE")

    Nt = length(KE.times)
    case[:times]      = KE.times
    case[:iterations] = collect(0:Nt-1)
    case[:KE]         = [KE[i][1, 1, 1]  for i in 1:Nt]
    case[:APE]        = [APE[i][1, 1, 1] for i in 1:Nt]
    case[:E]          = case[:KE] .+ case[:APE]
    case[:μ₀]         = μ₀

    λ, λ_err = seiche_damping_rate(case[:E])
    case[:λ]     = λ
    case[:λ_err] = λ_err

    ω = seiche_frequency(case[:KE], case[:times])
    case[:ω] = ω

    param = internal_seiche_parameters()
    c₀ = seiche_vertical_mode(0, param).cq
    c₁ = seiche_vertical_mode(param.q, param).cq
    Δt = length(case[:times]) > 1 ? case[:times][2] - case[:times][1] : NaN

    # ω = k_eff c₁ for the excited mode, so k_eff = ω / c₁ and μ₀_eff = k_eff c₀ Δt.
    case[:μ₀_eff] = isfinite(ω) ? (ω / c₁) * c₀ * Δt : NaN

    return case
end

"""
    seiche_damping_rate(E; skip=0.05, use=0.6)

Return the damping rate per time step of a decaying mode, and its standard error.

$(SIGNATURES)

# Arguments
- `E`: Total perturbation energy sampled once per time step
- `skip`: Fraction of the record to discard at the start, while the barotropic mode is being filtered out
- `use`: Fraction of the record to fit over

# Returns
- `(λ, λ_err)`: the amplitude damping rate per step and its standard error

The fit is a least squares line through `log E` against iteration number, halved because the energy is
quadratic in the amplitude. The opening fraction is discarded because the initial condition projects a
little onto the barotropic mode, which the averaging kernel removes within a handful of steps.
"""
function seiche_damping_rate(E; skip = 0.05, use = 0.6)
    N  = length(E)
    i₀ = max(2, round(Int, skip * N))
    i₁ = min(N, i₀ + round(Int, use * N))

    y = @view E[i₀:i₁]
    any(e -> !isfinite(e) || e ≤ 0, y) && return (NaN, NaN)

    x = collect(Float64, (i₀:i₁) .- 1)
    ly = log.(y)
    n  = length(x)
    x̄, ȳ = sum(x)/n, sum(ly)/n
    Sxx = sum((x .- x̄).^2)
    slope = sum((x .- x̄) .* (ly .- ȳ)) / Sxx

    residual = ly .- (ȳ .+ slope .* (x .- x̄))
    s² = sum(residual.^2) / (n - 2)

    return (slope / 2, sqrt(s² / Sxx) / 2)
end

"""
    seiche_frequency(KE, times)

Return the angular frequency of the mode from the oscillation of its kinetic energy.

$(SIGNATURES)

# Arguments
- `KE`: Kinetic energy time series
- `times`: Corresponding times

# Returns
- Angular frequency [rad/s], or `NaN` if too few oscillations are present to measure one

Kinetic and available potential energy exchange twice per wave period, so the kinetic energy oscillates
at `2ω` about its decaying mean. The frequency is taken from the mean spacing of its minima, which is
robust to the decay without needing to detrend.
"""
function seiche_frequency(KE, times)
    length(KE) < 8 && return NaN

    # Minima of the kinetic energy, one every half period of the KE oscillation, i.e. every quarter
    # period of the wave.
    minima = Int[]
    for i in 2:length(KE)-1
        KE[i] < KE[i-1] && KE[i] < KE[i+1] && push!(minima, i)
    end
    length(minima) < 2 && return NaN

    Δ = [times[minima[i+1]] - times[minima[i]] for i in 1:length(minima)-1]
    quarter_period = sum(Δ) / length(Δ)

    return 2π / (4 * quarter_period)
end

"""
    internal_seiche_sweep(folder, timestepper, free_surface, μ₀s)

Load a sweep over `μ₀` and return the measured damping rates.

$(SIGNATURES)

# Returns
- Named tuple `(; μ₀, μ₀_eff, λ, λ_err)`, each a vector over the runs that could be loaded

Runs whose output is missing are skipped rather than raising, so a partially completed sweep still
plots.
"""
function internal_seiche_sweep(folder, timestepper, free_surface, μ₀s)
    μ₀, μ₀_eff, λ, λ_err = Float64[], Float64[], Float64[], Float64[]

    for μ in μ₀s
        case = try
            load_internal_seiche(folder, timestepper, free_surface, μ)
        catch err
            @warn "skipping μ₀ = $μ for $free_surface" exception = err
            continue
        end
        push!(μ₀, case[:μ₀]); push!(μ₀_eff, case[:μ₀_eff])
        push!(λ, case[:λ]);   push!(λ_err, case[:λ_err])
    end

    return (; μ₀, μ₀_eff, λ, λ_err)
end

"""
    internal_seiche_order(folder, timestepper, free_surface, μ₀s)

Measure the convergence order of a sweep by comparing against the exact modal solution.

$(SIGNATURES)

# Returns
- Named tuple `(; μ₀, error, order)`, `order` being the successive slopes of `log error` against `log Δt`

The exact solution of the linearised modal system decays not at all, so the error after a fixed
*physical* time is dominated by the scheme's own amplitude and phase error. Comparing across a `μ₀`
sweep at fixed end time therefore returns the global order: one for the forward--backward substep, two
for the three-stage one.
"""
function internal_seiche_order(folder, timestepper, free_surface, μ₀s)
    param = internal_seiche_parameters()
    μ₀, err = Float64[], Float64[]

    for μ in μ₀s
        case = try
            load_internal_seiche(folder, timestepper, free_surface, μ)
        catch
            continue
        end
        # Energy is conserved exactly by the continuous problem, so its drift is the scheme's error.
        E = case[:E]
        push!(μ₀, μ)
        push!(err, abs(E[end] / E[1] - 1))
    end

    order = [log(err[i]/err[i+1]) / log(μ₀[i]/μ₀[i+1]) for i in 1:length(μ₀)-1]

    return (; μ₀, error = err, order)
end
