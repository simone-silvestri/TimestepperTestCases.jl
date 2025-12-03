```@meta
CurrentModule = TimestepperTestCases
```

# TimestepperTestCases.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://simone-silvestri.github.io/TimestepperTestCases.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://simone-silvestri.github.io/TimestepperTestCases.jl/dev/)
[![Build Status](https://github.com/simone-silvestri/TimestepperTestCases.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/simone-silvestri/TimestepperTestCases.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/simone-silvestri/TimestepperTestCases.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/simone-silvestri/TimestepperTestCases.jl)

**TimestepperTestCases.jl** is a Julia package that provides test cases and diagnostics for evaluating time-stepping schemes in free-surface ocean models. This repository complements the paper "A Low-Storage Runge-Kutta Framework for Non-Linear Free-Surface Ocean Models" by providing:

- Three idealized test cases spanning linear to fully turbulent regimes
- Exact variance budget diagnostics for quantifying numerical mixing
- Linear stability analysis tools
- Notebooks for reproducing figures and analysis from the paper

## Overview

This package implements test cases designed to highlight differences in performance between time discretization schemes, particularly comparing low-storage Runge-Kutta (RK) methods with Adams-Bashforth (AB2) schemes. The test cases include:

1. **Internal Tide**: A linear test case isolating temporal discretization effects
2. **Idealized Coastal Baroclinic Adjustment**: A nonlinear case demonstrating submesoscale variability
3. **Equilibrated Channel Flow**: A long-term equilibrated configuration showing interaction between numerical and physical mixing

## Key Features

- **Exact Variance Budgets**: Diagnostics that are exact locally in both space and time, enabling quantification of numerical mixing from time discretization
- **Multiple Timesteppers**: Support for `QuasiAdamsBashforth2`, `SplitRungeKutta3`, and `SplitRungeKutta4`
- **Free Surface Options**: Both implicit and split-explicit free surface formulations
- **Comprehensive Diagnostics**: Reference potential energy (RPE), available potential energy (APE), kinetic energy, and numerical diffusivity diagnostics

## Quick Links

- [Getting Started](@ref) - Installation and basic usage
- [Experiments](@ref) - Detailed descriptions of the three test cases
- [Diagnostics](@ref) - Understanding variance budgets and numerical mixing
- [Stability Analysis](@ref) - Linear stability analysis tools
- [Notebooks](@ref) - Reproducing figures and analysis
- [API Reference](@ref) - Complete function documentation

## Citation

If you use this package in your research, please cite the associated paper:

```bibtex
@article{silvestri2025timestepper,
  title={A Low-Storage Runge-Kutta Framework for Non-Linear Free-Surface Ocean Models},
  author={Silvestri, Simone and others},
  journal={Journal of Advances in Modeling Earth Systems},
  year={2025}
}
```

## Related Packages

This package builds on [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl), a fast, friendly, flexible ocean-flavored fluid dynamics package for CPUs and GPUs.
