# Time-discretization theory

`TimestepperTheory` is the Julia package in `lib/TimestepperTheory` for the amplification factors, barotropic averaging filters, and coupled barotropic–baroclinic analysis in Sections 2–4 of the manuscript. It runs independently of the ocean simulations.

The [example guide](theory_example.md) explains how to run and interpret the experiments. The [manuscript example](generated/time_discretization.md) shows the computations and all nine figures, using the manuscript terminology, legends, and colors. Plotting commands are hidden in the rendered example; the complete executable source is `lib/TimestepperTheory/examples/time_discretization.jl`.

## Modal state and normalization

The modal state is `[u₀, …, uQ, h₀, …, hQ]`, where velocity is normalized by √(gH) and total pressure by gH. The barotropic Courant number is μ₀ = k√(gH)Δt, and the dimensionless stratification is ε = N²H/g. The barotropic mode is included in addition to the requested number of baroclinic modes.

```@example theory
using TimestepperTheory

modes = vertical_modes(0.05, 2)
configuration = scheme_configuration("WRK3-SE")
filter = averaging_filter(configuration.averaging)
M = evolution_matrix(2.0, modes; configuration.scheme, configuration.forcing, configuration.solver, filter)
size(M)
```

The modal fields are `x`, `A`, `α`, `ψ`, `χ`, and `Π`; stage matrices are returned as `A`, `B`, and `M`. Filter outputs use `Δτ` and `w`, and barotropic outputs use `W`, `W̃`, `ūₚ`, and `ūₚ⁺`. The SRK3 blended matrices use χᵀχ ≈ 1.

## Computations

The example covers outer amplification factors, barotropic damping and averaging weights, stage and evolution matrices, direct and indirect coupling corrections, slow-forcing reconstruction, stability, and global temporal accuracy. Its final precision check distinguishes small growth from round-off.

Use the [theory API](theory_api.md) for the arguments, normalization, and returned quantities. The [simulation experiments](experiments/internal_tide.md) and [variance-budget diagnostics](diagnostics.md) describe the Oceananigans test cases.

## Reproduce the figures

From the repository root:

```sh
julia --project=lib/TimestepperTheory/examples -e 'using Pkg; Pkg.develop(path="lib/TimestepperTheory"); Pkg.instantiate()'
julia --project=lib/TimestepperTheory/examples lib/TimestepperTheory/examples/time_discretization.jl
```

PDF, PNG, and CSV files are written to `lib/TimestepperTheory/output`. The default run uses the manuscript grids; `--quick` selects coarse grids for an execution check.
