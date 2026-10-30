# TimestepperTheory.jl

Julia implementation of the time-discretization theory in Sections 2–4 of *A conservative Runge–Kutta framework for nonlinear free-surface ocean models*. The package constructs the modal evolution matrices and reproduces the manuscript's amplification, filtering, coupling, stability, and temporal-accuracy experiments. It runs independently of the ocean simulations and requires no simulation data or GPU.

## Read the example

The [Literate example](examples/time_discretization.jl) follows the analysis from the outer stability polynomials to the coupled barotropic–baroclinic system. Equations introduce each calculation, and the numerical experiments appear directly in the script. Figure layouts, colors, legends, and output are collected in [time_discretization_figures.jl](examples/time_discretization_figures.jl).

The [example guide](../../docs/src/theory_example.md) describes the parameters and outputs. The rendered example is part of the repository's Documenter site.

## Run and render

Run these commands from the monorepo root with Julia 1.12:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop([PackageSpec(path=pwd()), PackageSpec(path="lib/TimestepperTheory")]); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Open `docs/build/generated/time_discretization.html`. The build executes the full example and renders its computations, equations, and nine figures. PDF, PNG, and CSV outputs are saved in `lib/TimestepperTheory/output`.

To run only the computations with the same environment:

```sh
julia --project=docs lib/TimestepperTheory/examples/time_discretization.jl
```

Append `--quick` for a coarse execution check. Quick and full runs write the same filenames; use a full run for manuscript comparisons. To render documentation using existing figures, run `julia --project=docs docs/make.jl --reuse-figures`.

## Find the mathematics

```text
src/
├── algorithms.jl             # Schemes, solvers, filters, and forcing reconstructions
├── amplification_factors.jl  # Outer stability polynomials and the physical AB2 root
├── vertical_modes.jl         # x, A, α, ψ, χ, Π, and exact modal evolution
├── averaging_filters.jl      # Averaging weights and centered moments
├── barotropic_response.jl    # W, filtered transport, and polynomial forcing outputs
├── stage_matrices.jl         # Aₘ, Bₘ, Mₘ, and the complete timestep
└── diagnostics.jl            # Uncoupled evolution, V, eigenvalue corrections, and global error
```

The modal state is `y = [u₀, …, uQ, h₀, …, hQ]`, with velocity divided by √(gH) and total pressure by gH. The parameter `ε = N²H/g` measures stratification, and `μ₀ = k√(gH)Δt` is the barotropic Courant number. The requested count Q excludes the barotropic mode.

Numerical choices are concrete Julia objects: `WickerSkamarock3()`, `ShuOsher3()`, `ModifiedRungeKutta4(a=23 // 100)`, and `AdamsBashforth2()` select the outer scheme. Barotropic solvers are `ForwardBackward()`, `RungeKutta3()`, `ExactSolver()`, and `ImplicitSolver()`. Averaging uses `NoFilter()`, `F1Filter()`, `F2Filter()`, or `F3Filter()`.

Forcing reconstruction uses `FrozenForcing()`, `SSPForcing()`, `FinalStageQuadraticForcing()`, or `ProgressiveForcing()`. The time-dependent reconstructions accept `sampling=MidpointSampling()` by default; `ContinuousSampling()` requires an exact, unfiltered barotropic solve. The solver's polynomial outputs use local stage time s = t′/stage_duration; reconstruction coefficients include the corresponding powers of the stage duration.

`stage_matrices` returns `A[m] = Aₘ`, `B[m] = Bₘ`, and `M[m] = Mₘ`, with M₀ = I used internally. The modal weights retain the finite-basis truncation error in χᵀχ, while the SRK3 blended matrices use χᵀχ ≈ 1. AB2 uses an augmented state containing the previous timestep, initialized exactly in the accuracy experiment.

## Reproduce the comparisons

The example starts with outer amplification factors, barotropic damping, and averaging weights. It then compares full-matrix growth with direct and indirect eigenvalue corrections, tests forcing reconstruction, scans stability for the six manuscript configurations, and measures global temporal error.

The filtering experiments use 64 nominal substeps; the coupled experiments use 48. Coupling decompositions use 128-bit arithmetic, and forcing and accuracy experiments use 256 bits. Stability scans begin in Float64 and recompute growth below 10⁻¹² in magnitude at 128 bits. A separate `precision.csv` records the weak-stratification check at both precisions.

Two parameter conventions need attention when comparing with the manuscript captions: the forcing-stratification experiment uses μ₀ = 2, while the local caption gives 0.1; the accuracy experiment uses k√(gH)T = 64, while the caption writes kc₀T = 64. Set `μ₀` in the forcing-reconstruction experiment to reproduce the first caption, or pass `horizon=64/modes.α[1]` to `global_error` for the second interpretation. The manuscript files are separate from this package.

## Test

```sh
julia --project=lib/TimestepperTheory -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

The tests check amplification factors, filter moments, modal identities, analytical coupling corrections, temporal convergence, and matrix evolution against independent stage-by-stage updates. Near an eigenvalue collision, the perturbation expansion can become singular; the full-matrix comparison and eigenvector overlap help interpret those samples.
