# TimestepperTheory.jl

Time discretization and linear barotropic–baroclinic coupling analysis for Sections 2–4 of the timestepper manuscript. This package lives in the TimestepperTestCases.jl monorepo and runs independently of the ocean simulations.

Read the [example guide](../../docs/src/theory_example.md) for execution, parameter choices, and output files. The repository documentation also includes the rendered example and the complete theory API.

## Run

From the monorepo root:

```sh
julia --project=lib/TimestepperTheory -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=lib/TimestepperTheory/examples -e 'using Pkg; Pkg.develop(path="lib/TimestepperTheory"); Pkg.instantiate()'
julia --project=lib/TimestepperTheory/examples lib/TimestepperTheory/examples/time_discretization.jl
```

The Literate example `examples/time_discretization.jl` walks through the outer amplification factors, filter weights and damping, modal basis, stage matrices, and coupled experiments of Sections 2–4. The computations appear directly in the example; plotting commands use Literate’s `#src` marker to stay out of the rendered page. The shared file-output helper is in `examples/time_discretization_figures.jl`. Both computation and plotting are written in Julia.

To render the example with its generated figures:

```sh
julia --project=lib/TimestepperTheory/examples lib/TimestepperTheory/examples/render_time_discretization.jl
```

Open `lib/TimestepperTheory/output/time_discretization.html`. The reproduction script saves PDF, PNG, and CSV files in `lib/TimestepperTheory/output`. Add `--quick` for a coarse execution check with nine Courant-number samples; manuscript figures use the full grids. The numerical core needs LinearAlgebra, GenericSchur, and DocStringExtensions; CairoMakie, LaTeXStrings, and Literate belong to the separate example environment. Julia 1.10 or newer is required.

## Algorithm choices

Numerical choices use concrete objects and multiple dispatch. For example, `averaging_filter(F3Filter())` builds the F3 weights,
and `evolution_matrix(μ₀, modes; scheme=ShuOsher3(), forcing=SSPForcing(), solver=RungeKutta3())` constructs the SRK3 matrix.

The outer schemes are `WickerSkamarock3()`, `ShuOsher3()`, `ModifiedRungeKutta4()`, and `AdamsBashforth2()`. Barotropic solvers
are `ForwardBackward()`, `RungeKutta3()`, `ExactSolver()`, and `ImplicitSolver()`. Filters are `NoFilter()`, `F1Filter()`,
`F2Filter()`, and `F3Filter()`. Forcing reconstruction uses `FrozenForcing()`, `FinalStageQuadraticForcing()`,
`ProgressiveForcing()`, or `SSPForcing()`. The manuscript configuration helper returns these objects in its `scheme`, `solver`,
`forcing`, and `averaging` fields.

## Coordinates and time

The state is `y = [u₀, …, uQ, h₀, …, hQ]`. Velocity is divided by √(gH), and total pressure by gH. `μ₀ = k√(gH)Δt` is the barotropic Courant number and `ε = N²H/g` is the dimensionless stratification. Mode zero is included in addition to the requested `Q` baroclinic modes. Modal fields use `x`, `A`, `α`, `ψ`, `χ`, and `Π`; the quasi-Adams–Bashforth and MRK4 parameters are `ξ` and `a`.

The barotropic response uses local stage time s = t′/stage_duration. Its three forcing outputs correspond to 1, s, and s². Reconstruction weights include the stage-duration powers, so multiplying the outputs by reconstructed coefficients gives the same physical polynomial as a reconstruction in t′. `sampling=MidpointSampling()` holds the polynomial value fixed during a substep. `sampling=ContinuousSampling()` integrates the polynomial continuously and requires `solver=ExactSolver()` with an unfiltered barotropic solve.

The implicit configuration uses the pressure predictor and implicit surface solve at each WRK3 stage. The AB2 matrix includes the previous timestep in an augmented state; the accuracy calculation initializes that history from the exact solution.

## Mathematical implementation

| File | Quantities |
| --- | --- |
| `src/vertical_modes.jl` | Roots x tan(x)=ε; normalization; α, ψ, χ, Π; exact modal evolution |
| `src/averaging_filters.jl` | F1/F2/F3 averaging weights and moments |
| `src/barotropic_response.jl` | FB, RK3, exact and implicit barotropic evolution matrices; velocity and continuity transport |
| `src/stage_matrices.jl` | Aₘ, Bₘ, Mₘ; frozen, quadratic, progressive and SSP-weighted forcing; complete RK and AB2 matrices |
| `src/diagnostics.jl` | Uncoupled reference, coupling coefficients Vⱼₗ, direct/indirect corrections, frozen predictions, global error |
| `test/direct_steps.jl` | Independent state-by-state implementation of reconciliation and transport updates |

`stage_matrices` returns `A[m] = Aₘ`, `B[m] = Bₘ`, and `M[m] = Mₘ`; the initial matrix `M₀ = I` is used internally. Frozen forcing enters `B`; polynomial reconstruction contributes to `A` through the preceding evolution matrices. The final SRK3 stage uses the β-weighted forcing, followed by the Shu–Osher blend. These blended matrices use the approximation χᵀχ ≈ 1.

The finite basis retains the exact modal weights: χᵀχ is slightly less than one. The code does not renormalize χ. Consequently, I−χχᵀ is the truncated depth-average subtraction operator, not an exact orthogonal projector. The tests measure the omitted weight and verify its decay with mode count.

## Manuscript experiments

| Output | Experiment |
| --- | --- |
| `amplification_factor` | Physical AB2 root and RK2/RK3/RK4/MRK4 amplification magnitudes from Section 2 |
| `barotropic-filter-damping` | Exact, FB, and RK3 substeps with F1/F2/F3 and the implicit reference; 64 nominal substeps |
| `averaging-weights` | Output and continuity-transport weights, with centered moments M₂, M₃, M₄; 64 nominal substeps |
| `coupling_decomposition` | Full growth, complete second-order perturbation, retained direct/barotropic indirect sum, and leading prediction; ε=0.0005/0.05, Q=8 |
| `damping_balance` | Direct, indirect and total growth for exact unfiltered, exact/F3, RK3/F3, FB/F3 and FB/F1 barotropic solves; ε=0.01, Q=8 |
| `forcing_stratification` | WRK3 frozen/quadratic/progressive forcing at μ₀=2, Q=4, ε=0.64/4ᵏ; CSV also includes finite-substep comparisons |
| `mrk4_forcing` | MRK4 direct and full-matrix magnitude corrections divided by θ₁⁴; a=0.23, ε=0.05, Q=2, RK3/F3 |
| `schemes` | Full-spectrum growth for all six configurations at ε=0.05/0.15 and Q=8/32 |
| `order_schemes` | Fixed-horizon global error at dimensionless time 64, ε=0.05, Q=2, with exact AB2 startup |
| `precision.csv` | Repeat weak-stratification FB/F3 growth at 128 and 256 bits for Q=8/32 |

WRK3-UP shares the linear WRK3-SE matrix because the tracer reconstruction does not enter this linearized problem. WRK3-SE uses final-stage quadratic forcing; SRK3-SE uses the final SSP-weighted forcing; MRK4-SE uses progressive forcing with a least-squares quadratic at its fourth stage. F1 cases use FB; other split-explicit cases use RK3/F3. The Section 3 filter figures use 64 nominal substeps; coupled experiments use 48.

Forcing and accuracy experiments use 256-bit arithmetic; coupling decompositions use 128 bits. Following the manuscript calculation, full-spectrum scans use Float64 and recompute growth below 10⁻¹² in magnitude at 128 bits. Each CSV row records the precision used. Only the residual matrix in the global-error calculation is converted to Float64 for its spectral norm. Figure rendering converts final coordinates to Float64. The norm is the operator 2-norm in the modal variables (u, h/α).

The perturbation expansion is singular at an exact eigenvalue collision. The CSV records the overlap used to identify the selected full-matrix eigenmode; a poorly matched or nearly resonant sample needs interpretation. Finite scans establish sampled linear stability, not nonlinear stability or a proof over a continuous interval. Continuously integrated and midpoint-held forcing are separate experiments.

## Manuscript settings checked on 29 October 2026

The forcing-stratification plot uses μ₀=2 as requested in the figure discussion; the local manuscript caption still says 0.1. Set `reference_phase` in the forcing-reconstruction experiment to reproduce that caption. The accuracy experiment follows the archived calculation at k√(gH)T=64. The local caption instead writes kc₀T=64, which differs if c₀ denotes the exact barotropic modal speed. To use that interpretation, pass `horizon=64/modes.α[1]` to `global_error`. Neither manuscript source nor figure assets are modified by this package.

## Verification

The standalone package includes analytical checks for amplification factors and filter weights, alongside the modal and coupled-stage tests. The MRK4 direct and full-matrix magnitude corrections agree with the archived 63-point curves to within 1.1×10⁻¹⁷ after scaling by θ₁⁴. These comparisons use separate direct stage implementations in the test suite and archived manuscript calculations for numerical cross-checks.
