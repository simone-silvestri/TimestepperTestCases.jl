# Time-discretization theory

`TimestepperTheory` is the Julia package in `lib/TimestepperTheory` for the amplification factors, barotropic averaging
filters, and coupled barotropic–baroclinic analysis in Sections 2–4 of the manuscript. It runs on a CPU, independently of
Oceananigans and of the simulation outputs.

The [manuscript example](generated/time_discretization.md) shows the computations and all nine figures, with the manuscript
terminology, legends, and colors. It is a Literate.jl script: prose explains each experiment, the displayed code performs
the computations, and lines marked `#src` (plotting and file output) are omitted from the rendered page. The source is
`lib/TimestepperTheory/examples/time_discretization.jl`, its plotting helpers are in `time_discretization_figures.jl`, and
the numerical implementation is in `lib/TimestepperTheory/src/`. The [theory API](theory_api.md) documents the arguments,
normalization, and returned quantities.

## Modal state and normalization

The modal state is `[u₀, …, uQ, h₀, …, hQ]`, where velocity is normalized by √(gH) and total pressure by gH. The barotropic
Courant number is μ₀ = k√(gH)Δt, and the dimensionless stratification is ε = N²H/g. The modal count Q excludes the
barotropic mode, which is always included in addition to the requested baroclinic modes.

```@example theory
using TimestepperTheory

modes = vertical_modes(0.05, 2)
configuration = scheme_configuration("WRK3-SE")
filter = averaging_filter(configuration.averaging)
M = evolution_matrix(2.0, modes; configuration.scheme, configuration.forcing, configuration.solver, filter)
size(M)
```

The modal fields are `x`, `A`, `α`, `ψ`, `χ`, and `Π`; stage matrices are returned as `A`, `B`, and `M`. Filter outputs use
`Δτ` and `w`, and barotropic outputs use `W`, `W̃`, `ūₚ`, and `ūₚ⁺`. The SRK3 blended matrices use χᵀχ ≈ 1.

## Run the example

From the repository root, with Julia 1.12:

```sh
julia --project=lib/TimestepperTheory/examples -e 'using Pkg; Pkg.develop(path="lib/TimestepperTheory"); Pkg.instantiate()'
julia --project=lib/TimestepperTheory/examples lib/TimestepperTheory/examples/time_discretization.jl
```

No simulation data or GPU is required. Appending `--quick` reduces the sampled Courant-number grids of the coupled
experiments for an execution check; the amplification, filter, and stratification calculations keep their grids. Quick
and full runs write the same filenames, so run the full example again before using its figures in the manuscript. The
documentation build (`docs/make.jl`) executes the full example unless `--reuse-figures` is passed, in which case it uses
the existing PNG files and fails if one is missing.

## Experiments and outputs

The example covers outer amplification factors, barotropic damping and averaging weights, stage and evolution matrices,
direct and indirect coupling corrections, slow-forcing reconstruction, stability, and global temporal accuracy. Each
figure is saved as PDF and PNG, with a CSV of its numerical values, in `lib/TimestepperTheory/output`.

| Output stem | Computation |
| --- | --- |
| `amplification_factor` | Magnitude departure of the physical AB2 amplification factor and the RK stability polynomials for unit-speed advection |
| `barotropic-filter-damping` | Barotropic damping for exact, forward–backward, and RK3 substeps with F1, F2, and F3, compared with the implicit solver |
| `averaging-weights` | Averaging and continuity-transport weights, and their centered moments |
| `coupling_decomposition` | Full-matrix growth compared with direct and indirect eigenvalue corrections and the leading analytical prediction |
| `damping_balance` | Direct correction, barotropic indirect correction, and total growth for the barotropic solvers and filters |
| `forcing_stratification` | Dependence of coupling corrections on stratification for frozen and progressive forcing |
| `mrk4_forcing` | MRK4 direct and full-matrix magnitude corrections normalized by θ₁⁴ |
| `schemes` | Maximum growth over the full spectrum for the six manuscript configurations |
| `order_schemes` | Global temporal error at fixed final time |

`precision.csv` records repeated weak-stratification calculations at 128 and 256 bits; its final check distinguishes small
growth from round-off. It has no separate figure.

## Parameters

The filtering experiments use 64 nominal substeps; the coupled experiments use 48. Coupling decompositions use 128-bit
arithmetic, while forcing and temporal-accuracy experiments use 256 bits. The full-spectrum stability scan starts in Float64
and recomputes growth smaller than 10⁻¹² in magnitude at 128 bits.

The forcing-stratification experiment uses μ₀ = 2 and four baroclinic modes. Its CSV includes continuously integrated
forcing and midpoint-held forcing, while its plotted comparison uses continuous forcing. The MRK4 experiment uses a = 0.23,
ε = 0.05, two baroclinic modes, and RK3 substeps with F3.

The global-error experiment reaches the dimensionless time k√(gH)T = 64. It uses the operator norm in energy-scaled
coordinates [u; h/α], and starts AB2 with exact previous-step history. This final time uses the reference speed √(gH)
rather than the exact barotropic modal speed c₀.

## Reading the coupling results

The uncoupled evolution matrix advances each vertical mode independently. The difference between the coupled and
uncoupled matrices gives the coupling coefficients in the normalized left and right eigenvector basis. The direct
correction acts within the selected mode and propagation direction; the indirect correction includes transfers through
other modes.

The returned `indirect` correction includes the two barotropic propagation directions, while `other` contains the remaining
indirect terms. The example compares their linearized magnitude corrections with the full-matrix growth. The `overlap` field
measures the eigenvector match used to identify the selected full-matrix eigenvalue. Near an eigenvalue collision, the
perturbation estimate requires care; sampled stability curves are not proofs of stability for every Courant number.

## Change an experiment

Edit the parameters in the corresponding block of `time_discretization.jl`, then rerun the example. Scheme, solver, filter,
and forcing choices are concrete Julia objects such as `WickerSkamarock3()`, `RungeKutta3()`, `F3Filter()`, and
`ProgressiveForcing()`, documented in the [theory API](theory_api.md).
