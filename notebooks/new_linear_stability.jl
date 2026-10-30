# Standalone validation of the new notebook cells (no CairoMakie/Oceananigans needed).
# Ports: filters, weights_from_substeps, substep integrators, substep_matrices, SplitExplicitBaro
# compute_A/compute_B, compute_evolution_matrix, amplification. Checks against matrix_fig panel (c).
using LinearAlgebra
const FT = Float64

# ---- filters (from notebook cell 1) ----
function low_dispersion_filter(τ::T; p=2, q=4, r=0.18927) where T
    τ₀ = (p + 2) * (p + q + 2) / (p + 1) / (p + q + 1)
    return (τ / τ₀)^p * (1 - (τ / τ₀)^q) - r * (τ / τ₀)
end
no_filter(τ::T) where T = τ == 2 ? 1 : 0

# ---- weights (notebook cell 2) ----
@inline function weights_from_substeps(substeps, averaging_kernel)
    τᶠ = range(FT(0), FT(2), length = substeps+1)
    Δτ = τᶠ[2] - τᶠ[1]
    averaging_weights = map(averaging_kernel, τᶠ[2:end])
    M★ = something(findlast(>(0), averaging_weights), firstindex(averaging_weights))
    trimmed_weights = averaging_weights[1:M★]
    trimmed_weights = trimmed_weights ./ sum(trimmed_weights)
    barycenter = sum(trimmed_weights .* (1:M★)) * Δτ
    Δτ = Δτ / barycenter
    transport_weights = [sum(trimmed_weights[i:M★]) for i in 1:M★] .* Δτ
    return FT(Δτ), map(FT, tuple(trimmed_weights...)), map(FT, tuple(transport_weights...))
end

# ---- μ₂=0 coefficient by bisection (matches Oceananigans low_dispersion_coefficient) ----
function averaging_second_moment(NS, r)
    Δτ, w, _ = weights_from_substeps(NS, τ -> low_dispersion_filter(τ; r))
    sum(w[m]*(m*Δτ - 1)^2 for m in eachindex(w))
end
function mu2_zero_coefficient(NS; rl=0.18927, rh=0.9)
    fl = averaging_second_moment(NS, rl)
    fh = averaging_second_moment(NS, rh)
    fl*fh > 0 && return 0.285
    for _ in 1:80
        rm = (rl+rh)/2; fm = averaging_second_moment(NS, rm)
        (fl*fm ≤ 0) ? (rh = rm) : (rl = rm; fl = fm)
    end
    return (rl+rh)/2
end

# ---- barotropic substep integrators (one fast substep of the (u0, η) oscillator, baroclinic h1m frozen) ----
struct FBSubstep end
struct RK3Substep end
struct FBRK32Substep{FT}; β₁::FT; β₂::FT; β₃::FT; end
FBRK32Substep(; β₁=0.516, β₂=0.532, β₃=0.331) = FBRK32Substep(β₁, β₂, β₃)

@inline function barotropic_substep(::FBSubstep, u0, η, h1m, γB, k, g, H, ψ)
    u0n = u0 - γB * im * k * g * (η - ψ * h1m)
    ηn  = η  - γB * im * k * H * u0n
    return u0n, ηn
end
@inline function barotropic_substep(::RK3Substep, u0, η, h1m, γB, k, g, H, ψ)
    R(U, E) = (-im*k*g*(E - ψ*h1m), -im*k*H*U)
    a1,b1 = R(u0, η);   u1,e1 = u0 + γB/3*a1, η + γB/3*b1
    a2,b2 = R(u1, e1);  u2,e2 = u0 + γB/2*a2, η + γB/2*b2
    a3,b3 = R(u2, e2);  return u0 + γB*a3, η + γB*b3
end
@inline function barotropic_substep(ts::FBRK32Substep, u0, η, h1m, γB, k, g, H, ψ)
    β₁, β₂, β₃ = ts.β₁, ts.β₂, ts.β₃
    Ψ(U) = -im*k*H*U
    Φ(E) = -im*k*g*(E - ψ*h1m)
    e13 = η + γB/3*Ψ(u0);  es   = β₁*e13 + (1-β₁)*η;            u13 = u0 + γB/3*Φ(es)
    e12 = η + γB/2*Ψ(u13); ess  = β₂*e12 + (1-β₂)*η;           u12 = u0 + γB/2*Φ(ess)
    e1  = η + γB*Ψ(u12);   esss = β₃*e1 + (1-2β₃)*e12 + β₃*η
    return u0 + γB*Φ(esss), e1
end

# 2×2 homogeneous operator ℂ and forcing matrix D (column 2 = unit-h1m response) for a substep scheme
function substep_matrices(scheme, γB, k, g, H, ψ)
    ℂ = zeros(ComplexF64, 2, 2)
    D = zeros(ComplexF64, 2, 2)
    for (j, e) in enumerate(((1.0+0im, 0.0+0im), (0.0+0im, 1.0+0im)))
        u, η = barotropic_substep(scheme, e[1], e[2], 0.0+0im, γB, k, g, H, ψ)
        ℂ[1,j] = u; ℂ[2,j] = η
    end
    u, η = barotropic_substep(scheme, 0.0+0im, 0.0+0im, 1.0+0im, γB, k, g, H, ψ)
    D[1,2] = u; D[2,2] = η
    return ℂ, D
end

# ---- free-surface type carrying a barotropic substep scheme; w† (transport-corrected) coupling ----
struct SplitExplicitBaro{S}; substep::S; end

function compute_A(fs::SplitExplicitBaro, γ, k, g, N², H, ψ, ψ2=ψ; NS=NS, averaging_kernel=low_dispersion_filter)
    A = zeros(ComplexF64, 4, 4)
    Δτ, w, t = weights_from_substeps(NS, averaging_kernel)
    γB = γ * Δτ
    c² = N² * H^2 / π^2
    C  = γ * c² / g
    ℂ, D = substep_matrices(fs.substep, γB, k, g, H, ψ)
    M = length(w)
    ABT = sum([w[i]*ℂ^i for i in 1:M])
    ÃBT = sum([t[i]*ℂ^i for i in 1:M])
    A[1,1] = ABT[1,1]
    A[1,2] = ABT[1,2]
    A[1,4] = ABT[1,2] * ψ
    A[2,1] = ABT[2,1] - im*k*C*ψ*ψ2*ÃBT[1,1]
    A[2,2] = ABT[2,2] - im*k*C*ψ*ψ2*ÃBT[1,2]
    A[2,4] = ABT[2,2]*ψ - ψ - im*k*C*ψ^2*ψ2*ÃBT[1,2]
    A[3,3] = 1.0
    A[4,1] = im*k*C*ψ2*ÃBT[1,1]
    A[4,2] = im*k*C*ψ2*ÃBT[1,2]
    A[4,4] = 1.0 + im*k*C*ψ*ψ2*ÃBT[1,2]
    return A
end

function compute_B(fs::SplitExplicitBaro, γ, k, g, N², H, ψ, ψ2=ψ; NS=NS, averaging_kernel=low_dispersion_filter)
    B = zeros(ComplexF64, 4, 4)
    Δτ, w, t = weights_from_substeps(NS, averaging_kernel)
    γB = γ * Δτ
    c² = N² * H^2 / π^2
    C  = γ * c² / g
    ℂ, D = substep_matrices(fs.substep, γB, k, g, H, ψ)
    M = length(w)
    BBT = sum([w[i]*sum([ℂ^s*D for s in 0:i-1]) for i in 1:M])
    B̃BT = sum([t[i]*sum([ℂ^s*D for s in 0:i-1]) for i in 1:M])
    B[1,4] = BBT[1,2]
    B[2,3] = im * k * C * ψ
    B[2,4] = BBT[2,2] - im*k*C*ψ*ψ2*B̃BT[1,2]
    B[3,4] = - γ * g * im * k
    B[4,3] = -im * k * C
    B[4,4] = im * k * C * ψ2 * B̃BT[1,2]
    return B
end

# ---- RK3 baroclinic assembly (notebook cell 14, Val(3)) ----
function compute_evolution_matrix(free_surface, Δt, k, g, N², H; ψ = sqrt(2.0), ψ2 = ψ)
    γ₁ = Δt/3; A1 = compute_A(free_surface, γ₁, k, g, N², H, ψ, ψ2); B1 = compute_B(free_surface, γ₁, k, g, N², H, ψ, ψ2)
    γ₂ = Δt/2; A2 = compute_A(free_surface, γ₂, k, g, N², H, ψ, ψ2); B2 = compute_B(free_surface, γ₂, k, g, N², H, ψ, ψ2)
    γ₃ = Δt;   A3 = compute_A(free_surface, γ₃, k, g, N², H, ψ, ψ2); B3 = compute_B(free_surface, γ₃, k, g, N², H, ψ, ψ2)
    return A3 + B3*A2 + B3*B2*A1 + B3*B2*B1
end

# ---- parameters (cell 18) ----
g = 9.81; H = 6000.0; N² = 1e-4
Δx = 7812.0; k = 1/Δx; c0 = sqrt(g*H)
NS = 40
ψ = -sqrt(2); ψ2 = π^2/ψ/4
R_MU2 = mu2_zero_coefficient(NS)
mu2zero_filter(τ) = low_dispersion_filter(τ; r = R_MU2)

println("NS=$NS  μ₂=0 coefficient r=", round(R_MU2, digits=4),
        "  (μ₂ at r: ", round(averaging_second_moment(NS, R_MU2), sigdigits=2), ")")

function maxλ(fs, cfl; ker=low_dispersion_filter)
    Δt = cfl/(c0*k)
    # temporarily set the module-level averaging kernel via closure: use keyword
    M = compute_evolution_matrix_ker(fs, Δt, ker)
    maximum(abs.(eigvals(M)))
end
# helper that threads the kernel through
function compute_evolution_matrix_ker(fs, Δt, ker)
    γ₁ = Δt/3; A1 = compute_A(fs, γ₁, k, g, N², H, ψ, ψ2; averaging_kernel=ker); B1 = compute_B(fs, γ₁, k, g, N², H, ψ, ψ2; averaging_kernel=ker)
    γ₂ = Δt/2; A2 = compute_A(fs, γ₂, k, g, N², H, ψ, ψ2; averaging_kernel=ker); B2 = compute_B(fs, γ₂, k, g, N², H, ψ, ψ2; averaging_kernel=ker)
    γ₃ = Δt;   A3 = compute_A(fs, γ₃, k, g, N², H, ψ, ψ2; averaging_kernel=ker); B3 = compute_B(fs, γ₃, k, g, N², H, ψ, ψ2; averaging_kernel=ker)
    return A3 + B3*A2 + B3*B2*A1 + B3*B2*B1
end

schemes = (("FB + SM05",       SplitExplicitBaro(FBSubstep()),      low_dispersion_filter),
           ("RK3 + SM05",      SplitExplicitBaro(RK3Substep()),     low_dispersion_filter),
           ("FB-RK32 + SM05",  SplitExplicitBaro(FBRK32Substep()),  low_dispersion_filter),
           ("RK3 no-filter",   SplitExplicitBaro(RK3Substep()),     no_filter),
           ("FB-RK32 no-filter",SplitExplicitBaro(FBRK32Substep()), no_filter),
           ("FB μ₂=0",         SplitExplicitBaro(FBSubstep()),      mu2zero_filter),
           ("FB-RK32 μ₂=0",    SplitExplicitBaro(FBRK32Substep()),  mu2zero_filter))

println("\nmax|λ| at CFL = 2, 6, 10, 14:")
for (nm, fs, ker) in schemes
    vals = [maxλ(fs, c; ker) for c in (2.0, 6.0, 10.0, 14.0)]
    tag = maximum([maximum([maxλ(fs,c;ker) for c in range(0.2, 16, length=400)])]) # peak
    peak = maximum([maxλ(fs,c;ker) for c in range(0.2, 16, length=600)])
    println(rpad(nm,20), join([rpad(round(v,digits=4),9) for v in vals]), "  peak[0.2,16]=", round(peak,digits=4),
            peak ≤ 1.0005 ? "  STABLE" : "  UNSTABLE")
end

