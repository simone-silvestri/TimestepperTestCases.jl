# TimestepperTestCases.jl

This monorepo contains ocean-model test cases, buoyancy-variance diagnostics, and the time-discretization analysis accompanying *A conservative Runge–Kutta framework for nonlinear free-surface ocean models*.

## Time-discretization theory

`TimestepperTheory` implements the amplification factors, averaging filters, modal evolution matrices, and direct and indirect coupling corrections in Sections 2–4. It runs independently of the ocean simulations.

- [Theory overview](theory.md): modal state, normalization, and a first matrix calculation.
- [Example guide](theory_example.md): execution commands, parameters, and output files.
- [Manuscript example](generated/time_discretization.md): computations and all nine figures.
- [Theory API](theory_api.md): algorithm choices and numerical functions.

## Ocean simulations

`TimestepperTestCases` provides Oceananigans simulation setups and diagnostics for numerical mixing. The simulation and theory packages have separate Julia environments.

- [Getting started](getting_started.md)
- [Internal tide](experiments/internal_tide.md)
- [Re-entrant channel](experiments/channel.md)
- [Variance-budget diagnostics](diagnostics.md)
- [Simulation API](api_reference.md)
- [Data-driven Julia examples](data_examples.md)

## Build locally

From the repository root, using Julia 1.12:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop([PackageSpec(path=pwd()), PackageSpec(path="lib/TimestepperTheory")]); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The build executes the theory example and writes the site to `docs/build`. Pass `--reuse-figures` to use the existing figure outputs. The documentation workflow produces a downloadable site artifact; publishing is separate.
