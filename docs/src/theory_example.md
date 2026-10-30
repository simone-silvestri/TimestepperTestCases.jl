# Guide to the manuscript example

The [time-discretization example](generated/time_discretization.md) reproduces the amplification, filtering, and coupled modal analyses in Sections 2–4. It is a Julia script written with Literate.jl: prose explains each experiment, the displayed code performs the computations, and the figures show the resulting comparisons.

The source is `lib/TimestepperTheory/examples/time_discretization.jl`. The numerical implementation is in `lib/TimestepperTheory/src/`. Lines marked `#src` perform plotting and file output but are omitted from the rendered page. The shared `save_results` function writes each figure and its numerical data.

## Run and render

From the repository root, install the example dependencies once:

```sh
julia --project=lib/TimestepperTheory/examples -e 'using Pkg; Pkg.develop(path="lib/TimestepperTheory"); Pkg.instantiate()'
```

Execute the computations and render the standalone page:

```sh
julia --project=lib/TimestepperTheory/examples lib/TimestepperTheory/examples/time_discretization.jl
julia --project=lib/TimestepperTheory/examples lib/TimestepperTheory/examples/render_time_discretization.jl
```

Open `lib/TimestepperTheory/output/time_discretization.html`. The rendering command formats the source and links to the generated figures; it does not execute the computations. No simulation data or GPU is required. The first run also compiles Julia plotting dependencies.

For a coarse execution check, append `--quick` to the computation command. This reduces the sampled Courant-number grids in the coupled experiments; the amplification, filter, and stratification calculations retain their grids. Quick and full runs write the same filenames, so run the full example again before using its figures in the manuscript.

## Experiments and outputs

Each figure is saved as PDF and PNG, with a CSV containing its numerical values. All files are written to `lib/TimestepperTheory/output`.

| Output stem | Computation |
| --- | --- |
| `amplification_factor` | Magnitude departure of the physical AB2 amplification factor and the RK stability polynomials for unit-speed advection |
| `barotropic-filter-damping` | Barotropic damping for exact, forward–backward, and RK3 substeps with F1, F2, and F3, compared with the implicit solver |
| `averaging-weights` | Averaging and continuity-transport weights, and their centered moments |
| `coupling_decomposition` | Full-matrix growth compared with direct and indirect eigenvalue corrections and the leading analytical prediction |
| `damping_balance` | Direct correction, barotropic indirect correction, and total growth for the barotropic solvers and filters |
| `forcing_stratification` | Dependence of coupling corrections on stratification for frozen, final-stage quadratic, and progressive forcing |
| `mrk4_forcing` | MRK4 direct and full-matrix magnitude corrections normalized by θ₁⁴ |
| `schemes` | Maximum growth over the full spectrum for the six manuscript configurations |
| `order_schemes` | Global temporal error at fixed final time |

`precision.csv` records repeated weak-stratification calculations at 128 and 256 bits. It has no separate figure.

## Parameters and normalization

The [theory overview](theory.md) defines the modal state and its normalization. In the example, `phase` denotes the barotropic Courant number μ₀ = k√(gH)Δt, and `stratification` denotes ε = N²H/g. The modal count Q excludes the barotropic mode.

The filtering experiments use 64 nominal substeps; the coupled experiments use 48. Coupling decompositions use 128-bit arithmetic, while forcing and temporal-accuracy experiments use 256 bits. The full-spectrum stability scan starts in Float64 and recomputes growth smaller than 10⁻¹² in magnitude at 128 bits.

The forcing-stratification experiment uses μ₀ = 2 and four baroclinic modes. Its CSV includes continuously integrated forcing and midpoint-held forcing, while its plotted comparison uses continuous forcing. The MRK4 experiment uses a = 0.23, ε = 0.05, two baroclinic modes, and RK3 substeps with F3.

The global-error experiment reaches dimensionless time k√(gH)T = 64. It uses the operator norm in energy-scaled coordinates [u; h/α], and starts AB2 with exact previous-step history. This final time uses the reference speed √(gH), rather than the exact barotropic modal speed c₀.

## Reading the coupling results

The uncoupled evolution matrix advances each vertical mode independently. The difference between the coupled and uncoupled matrices gives the coupling coefficients in the normalized left and right eigenvector basis. The direct correction acts within the selected mode and propagation direction; the indirect correction includes transfers through other modes.

The returned `indirect` correction includes the two barotropic propagation directions, while `other` contains the remaining indirect terms. The example compares their linearized magnitude corrections with the full-matrix growth. The `overlap` field measures the eigenvector match used to identify the selected full-matrix eigenvalue. Near an eigenvalue collision, the perturbation estimate requires care; sampled stability curves are not proofs of stability for every Courant number.

## Change an experiment

Edit the parameters in the corresponding block of `time_discretization.jl`, then rerun the computation and rendering commands. Scheme, solver, filter, and forcing choices are concrete Julia objects such as `WickerSkamarock3()`, `RungeKutta3()`, `F3Filter()`, and `FinalStageQuadraticForcing()`. Their definitions and numerical outputs are documented in the [theory API](theory_api.md).

To build the complete repository documentation, use the commands in the repository README. `docs/make.jl` executes the full example by default; `--reuse-figures` uses the existing PNG files and fails if a required figure is missing.
