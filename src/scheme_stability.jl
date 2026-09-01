using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces: weights_from_substeps
using Oceananigans.TimeSteppers: SplitRungeKuttaTimeStepper, SSPRungeKuttaTimeStepper

#####
##### Choosing the baroclinic time step and the barotropic substep count from the linear theory
#####
##### The time step of every test case is set from the theoretical stability limit of the composition on the
##### first baroclinic mode, rather than from a cost match. Two Courant numbers are involved and both have to
##### be respected:
#####
##### - the baroclinic Courant number μ₁ = c₁ k Δt, which the outer Runge-Kutta scheme caps at the imaginary-
#####   axis limit θ★ of its amplification polynomial. This sets Δt.
##### - the barotropic substep Courant number c₀ k Δτ, which the sub-step integrator caps at its own limit.
#####   This sets the number of substeps, given Δt.
#####
##### `k` is the largest wavenumber the grid carries, π/Δx, since stability must hold for every resolved mode
##### and the grid-scale one is the binding constraint. The barotropic sub-cycle is the one place where the
##### distinction between that spectral wavenumber and the one a two-point difference actually delivers is
##### worth making -- see [`staggered_wavenumber`](@ref).
#####

"""
    amplification_polynomial(β, z)

Amplification factor of the restarted low-storage composition with coefficients `β`, that is the polynomial
`R(z)` such that `Ψⁿ⁺¹ = R(λΔt) Ψⁿ` for `G = λΨ`.

Each stage restarts from `Ψⁿ` and advances by `γₘ = Δt/βₘ`, so the polynomial is the nested product
`Rₘ = 1 + (z/βₘ) Rₘ₋₁`. For `β = (3, 2, 1)` this is the degree-three Taylor polynomial, and for
`β = (1/a, 3, 2, 1)` it is that polynomial plus `a z⁴/6`.
"""
@inline amplification_polynomial(β, z) = foldl((R, b) -> 1 + z * R / b, β; init = one(z))

"""
    imaginary_axis_limit(β; tolerance = 1e-6)

Largest `θ` for which `|R(iθ)| ≤ 1`, that is the stability limit of the composition on a purely oscillatory
mode. Returns `√3` for the three-stage scheme and `2√2` for the classical four-stage one.
"""
function imaginary_axis_limit(β; tolerance = 1e-6)
    amplifies(θ) = abs(amplification_polynomial(β, im * θ)) > 1

    lower, upper = 0.0, 0.1
    while !amplifies(upper) && upper < 100
        lower, upper = upper, 2upper
    end

    while upper - lower > tolerance
        middle = (lower + upper) / 2
        amplifies(middle) ? (upper = middle) : (lower = middle)
    end

    return lower
end

# The two-level scheme is not a one-step polynomial: its amplification factors are the two roots of the
# companion quadratic of Ψⁿ⁺¹ = Ψⁿ + Δt[(3/2 + χ) Gⁿ - (1/2 + χ) Gⁿ⁻¹], namely λ² - b λ + c = 0.
function quasi_ab2_limit(χ; tolerance = 1e-6)
    function amplifies(θ)
        z = im * θ
        b = 1 + (3/2 + χ) * z
        c = (1/2 + χ) * z
        discriminant = sqrt(b^2 - 4c)
        return max(abs((b + discriminant) / 2), abs((b - discriminant) / 2)) > 1 + 1e-12
    end

    lower, upper = 0.0, 0.1
    while !amplifies(upper) && upper < 100
        lower, upper = upper, 2upper
    end

    while upper - lower > tolerance
        middle = (lower + upper) / 2
        amplifies(middle) ? (upper = middle) : (lower = middle)
    end

    return lower
end

"""
    stability_limit(scheme)

Imaginary-axis stability limit `θ★` of a baroclinic time-stepping scheme, given either as the symbol used by
`HydrostaticFreeSurfaceModel` or as a materialized timestepper.
"""
stability_limit(β::Tuple) = imaginary_axis_limit(β)
stability_limit(ts::SplitRungeKuttaTimeStepper) = imaginary_axis_limit(ts.β)
stability_limit(::SSPRungeKuttaTimeStepper) = imaginary_axis_limit((3, 2, 1))  # same polynomial as WRK3

function stability_limit(scheme::Symbol)
    scheme === :QuasiAdamsBashforth2 && return quasi_ab2_limit(0.1)
    scheme === :SSPRungeKutta3       && return imaginary_axis_limit((3, 2, 1))
    scheme === :ModifiedRungeKutta4 && return imaginary_axis_limit(ModifiedRungeKutta4TimeStepper().β)

    name = string(scheme)
    startswith(name, "SplitRungeKutta") || throw(ArgumentError("Unknown scheme $scheme."))
    stages = parse(Int, name[length("SplitRungeKutta")+1:end])

    return imaginary_axis_limit(tuple(stages:-1:1...))
end

#####
##### Wave speeds
#####

"""
    first_baroclinic_speed(N², H)

Phase speed `c₁ = N H / π` of the first baroclinic mode of a constant stratification `N²` over a depth `H`.

This is the rigid-lid limit. The exact free-surface speed is `α₁ √(gH)` with `α₁` the root of `x tan x = ε`,
`ε = N²H/g`, which differs by `O(ε)`; the seiche case, whose whole point is that distinction, uses the exact
modal solver of [`seiche_vertical_mode`](@ref) instead.
"""
@inline first_baroclinic_speed(N², H) = sqrt(N²) * H / π

"""
    first_baroclinic_speed(N²::Function, H; levels = 2048)

WKB phase speed `c₁ = (1/π) ∫₋H⁰ N dz` for a depth-dependent stratification, with `N²` a function of `z`.

Reduces to `N H / π` for a constant `N`, and is the right estimate where the stratification is strongly
surface-intensified, as in the channel: taking the surface `N²` as if it were uniform would overstate `c₁` by
more than a factor of two and shorten the time step accordingly.
"""
function first_baroclinic_speed(N²::Function, H; levels = 2048)
    Δz = H / levels
    z  = range(-H + Δz/2, -Δz/2, length = levels)
    return sum(sqrt(max(N²(zᵢ), 0)) for zᵢ in z) * Δz / π
end

"""
    barotropic_speed(H; g = Oceananigans.defaults.gravitational_acceleration)

External gravity wave speed `c₀ = √(gH)`.
"""
@inline barotropic_speed(H; g = Oceananigans.defaults.gravitational_acceleration) = sqrt(g * H)

"""
    grid_wavenumber(Δx)

Largest wavenumber the grid carries, `π/Δx`. Stability must hold for every resolved mode, and the grid-scale
mode is the binding one, so this is the wavenumber that enters both Courant numbers.
"""
@inline grid_wavenumber(Δx) = π / Δx

"""
    staggered_wavenumber(Δx, Δy)

Largest wavenumber a *two-point staggered* difference delivers, `2√(Δx⁻² + Δy⁻²)`.

The symbol of the difference between two points a distance `Δx` apart is `2i sin(kΔx/2)/Δx`, which reaches
`2/Δx` at the grid-scale mode rather than the spectral `π/Δx` of [`grid_wavenumber`](@ref): a two-point
operator does not see the grid-scale mode at its own wavenumber. The barotropic sub-cycle advances `η` and
`(U, V)` with exactly such operators, so `2√(Δx⁻² + Δy⁻²)` is the spectral radius of the semi-discrete
barotropic problem divided by `c₀`, the square root over the two directions being the worst case over the
direction of propagation.

The two conventions agree to 11% on an isotropic grid, `2√2` against `π`, which is why the idealized cases are
indifferent to the choice; where the grid is strongly anisotropic they are not, and the spectral one is then
the more conservative of the two.
"""
@inline staggered_wavenumber(Δx, Δy) = 2 * sqrt(1 / Δx^2 + 1 / Δy^2)

#####
##### The two selections
#####

"""
    baroclinic_stability_limit(scheme; N², H, Δx)

The theoretical limit `Δt = θ★ / (c₁ k)` beyond which `scheme` amplifies the first baroclinic mode at the
grid-scale wavenumber `k = π/Δx`.
"""
function baroclinic_stability_limit(scheme; N², H, Δx)
    c₁ = first_baroclinic_speed(N², H)
    k  = grid_wavenumber(Δx)
    return stability_limit(scheme) / (c₁ * k)
end

"""
    advective_stability_limit(scheme; U, Δx, horizontal_dimensions = 2)

The limit `Δt = θ★ / (√d U k)` beyond which `scheme` amplifies advection at the grid scale.

Advection has the same structure as the wave problem: its semi-discrete eigenvalue is `i(u kₓ + v k_y)`, purely
imaginary for a centered operator, so the same `θ★` applies with the wave speed replaced by the advective one.
The `√d` accounts for direction: the eigenvalue is largest when the flow is diagonal to the grid, where
`|u|kₓ + |v|k_y = √2 U k` for a speed `U` in two dimensions. This is the same worst case that `k = π/Δx`
represents in wavenumber, and it agrees with the usual `|u|/Δx + |v|/Δy` sum.
"""
advective_stability_limit(scheme; U, Δx, horizontal_dimensions = 2) =
    stability_limit(scheme) / (sqrt(horizontal_dimensions) * U * grid_wavenumber(Δx))

"""
    baroclinic_timestep(scheme; N², H, Δx, U = nothing, horizontal_dimensions = 2,
                        safety = 0.7, reference_scheme = :SplitRungeKutta3)

Baroclinic time step of `scheme`: `safety` times the theoretical limit of the reference three-stage scheme,
scaled by the ratio of the imaginary-axis limits `θ★(scheme) / θ★(reference)`.

Anchoring on the reference scheme keeps every scheme on one operating point of one configuration, so that a
comparison between schemes is a comparison of the compositions and not of two independently chosen time steps.
The four-stage MRK4 composition earns `2.755/√3 = 1.591` over the three-stage schemes, and the two-level
scheme gives up `0.502/√3 = 0.29`.
"""
function baroclinic_timestep(scheme; N², H, Δx, U = nothing, horizontal_dimensions = 2,
                             safety = 0.7, reference_scheme = :SplitRungeKutta3)

    ratio = stability_limit(scheme) / stability_limit(reference_scheme)
    Δt_wave = safety * ratio * baroclinic_stability_limit(reference_scheme; N², H, Δx)

    isnothing(U) && return Δt_wave

    # Both limits carry the same factor θ★(scheme), so taking the smaller one still leaves every scheme at the
    # same fraction of its own limit -- it only changes which physics sets that limit.
    Δt_advective = safety * advective_stability_limit(scheme; U, Δx, horizontal_dimensions)

    return min(Δt_wave, Δt_advective)
end

"""
    baroclinic_timestep(scheme, reference_Δt; reference_scheme = :SplitRungeKutta3)

Time step of `scheme` scaled from a *measured* reference time step, for configurations where the theoretical
limit is not well posed -- a global domain, where `c₁` and `Δx` both vary and the binding `c₁ k` is a maximum
over the globe that needs the actual stratification field.

Only the ratio `θ★(scheme)/θ★(reference)` is used, which is the part of the criterion that transfers
regardless of the configuration.
"""
baroclinic_timestep(scheme, reference_Δt; reference_scheme = :SplitRungeKutta3) =
    reference_Δt * stability_limit(scheme) / stability_limit(reference_scheme)

"""
    barotropic_substeps(barotropic_scheme; H, Δx, Δt, averaging_kernel, safety = 0.7, granularity = 8,
                        wavenumber = grid_wavenumber(Δx))

Smallest number of barotropic substeps for which the substep Courant number `c₀ k Δτ` stays at or below
`safety` times the limit of the substep integrator, rounded up to a multiple of `granularity`.

`k` is `grid_wavenumber(Δx)` unless `wavenumber` is given, in which case `Δx` is not used. The near-global case
supplies [`staggered_wavenumber`](@ref) there, its grid being anisotropic enough for the two to disagree.

The substep size is not `Δt / substeps`: the averaging kernel spans a window wider than the baroclinic step,
and `weights_from_substeps` returns the fractional step size that goes with the requested count, so the
Courant number is evaluated on the step the sub-cycle actually takes. Rounding up to a multiple of 8 keeps
the count compatible with the kernels that require `substeps % 8 == 0`.

The limits are those of the substep integrators: `√3` for the three-stage Runge-Kutta substep and `1` for
forward-backward, the latter halved from its neutral `2` -- see [`substep_limit`](@ref).
"""
function barotropic_substeps(barotropic_scheme; H, Δx, Δt, averaging_kernel,
                             safety = 0.7, granularity = 8, maximum_substeps = 4096,
                             wavenumber = grid_wavenumber(Δx))

    c₀ = barotropic_speed(H)
    k  = wavenumber
    limit = safety * substep_limit(barotropic_scheme)

    substeps = granularity
    while substeps ≤ maximum_substeps
        fractional_Δt, _, _ = weights_from_substeps(Float64, substeps, averaging_kernel)
        c₀ * k * fractional_Δt * Δt ≤ limit && return substeps
        substeps += granularity
    end

    throw(ArgumentError("No substep count below $maximum_substeps keeps the barotropic Courant number " *
                        "under $limit for Δt = $Δt."))
end

"""
    substep_limit(barotropic_scheme)

Substep Courant number `c₀ k Δτ` at which the substep integrator is held, one per integrator.

For the three-stage Runge-Kutta substep this is its own limit, `√3`: the amplification is the truncated
exponential of the oscillator generator, it damps genuinely throughout the range (`|det| = 0.89` at `μ = 1.5`),
and `√3` is where that damping turns into growth.

For forward-backward the corresponding number would be `2`, but that limit describes neutral stability rather
than usable behaviour, so `1` is used instead. Forward-backward has unit determinant at *every* Courant
number -- it never damps, at any `μ` -- and at exactly `μ = 2` its matrix is defective (a double eigenvalue at
`-1`), so the state grows linearly in the substep count while `|λ| = 1`. The transient amplification is already
severe below the limit, with the hundred-step norm reaching 18 at `μ = 1.99`. Since all of the dissipation then
comes from the averaging filter rather than from the substep, the effective limit is not something this
analysis determines, and the conservative value is the honest one.
"""
@inline substep_limit(::ForwardBackwardScheme) = 1.0
@inline substep_limit(::RungeKutta3Scheme) = sqrt(3)

"""
    barotropic_courant(; H, Δx, Δt, substeps, averaging_kernel, wavenumber = grid_wavenumber(Δx))

Substep Courant number `c₀ k Δτ` that `substeps` actually delivers, for reporting alongside its limit. `k`
follows the same convention as in [`barotropic_substeps`](@ref).
"""
function barotropic_courant(; H, Δx, Δt, substeps, averaging_kernel, wavenumber = grid_wavenumber(Δx))
    fractional_Δt, _, _ = weights_from_substeps(Float64, substeps, averaging_kernel)
    return barotropic_speed(H) * wavenumber * fractional_Δt * Δt
end
