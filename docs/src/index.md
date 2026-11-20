# TimestepperTestCases.jl

This repository accompanies the manuscript *A conservative Runge–Kutta framework for nonlinear free-surface ocean models*.
It contains two Julia packages:

- `TimestepperTestCases` sets up the three Oceananigans experiments of the manuscript and diagnoses their numerical
  mixing from the buoyancy-variance budget.
- `TimestepperTheory`, in `lib/TimestepperTheory`, computes the amplification factors, barotropic averaging filters and
  coupled modal stability of Sections 2–4. It runs on a CPU without Oceananigans or simulation data.

Each experiment compares the same discretizations: Adams–Bashforth 2 (AB2) and the Runge–Kutta schemes WRK3, SRK3 and
MRK4, with split-explicit (SE) or implicit (IM) barotropic solvers. Start with [Getting started](getting_started.md).

## Experiments

### Internal tide

A two-dimensional M₂ tide over a Gaussian ridge in a stratified ocean, integrated for 40 days. Advection plays a minor
role, so the differences between runs come from the time discretization.
[Configuration and run](experiments/internal_tide.md) · [Manuscript figures](generated/internal_tide_figures.md)

![Final surface height, velocity and stratification of the internal tide](generated/figure-internal-tide.png)

### Re-entrant channel

A wind- and buoyancy-forced channel integrated to statistical equilibrium, where numerical mixing competes with the
parameterized vertical diffusivity.
[Configuration and run](experiments/channel.md) · [Manuscript figures](generated/channel.md)

![Numerical diffusivity of the channel](generated/channel_diffusivity.png)

### Quarter-degree global ocean

A five-year eORCA025 configuration on a GPU, used to measure numerical diffusivity, surface variability and cost.
[Configuration and run](experiments/global_ocean.md) · [Numerical diffusivity](generated/global_ocean_results.md) ·
[Surface velocity](generated/global_surface_results.md) · [Computational cost](generated/global_ocean.md)

![Surface speed and regional spectra of the global configuration](generated/global_surface_comparison.png)

## Theory, data and diagnostics

- [Time-discretization theory](theory.md) and its [manuscript example](generated/time_discretization.md), which draws
  Figures 1–9.
- [Reproducing the paper](data_examples.md): the Zenodo record of the simulation outputs and the script behind every
  manuscript figure.
- [Variance-budget diagnostics](diagnostics.md): how numerical diffusivity is diagnosed from the variance budget.
- API references for the [simulation package](api_reference.md) and the [theory package](theory_api.md).

## Build the documentation

From the repository root:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop([PackageSpec(path=pwd()), PackageSpec(path="lib/TimestepperTheory")]); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The build executes the theory example and writes the site to `docs/build`; `--reuse-figures` reuses the existing theory
figures. The data-driven pages are rendered without running their scripts and show the figures stored in `docs/figures`.
Successful builds of `main` are published to GitHub Pages.
