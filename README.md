# TimestepperTestCases.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://simone-silvestri.github.io/TimestepperTestCases.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://simone-silvestri.github.io/TimestepperTestCases.jl/dev/)
[![Build Status](https://github.com/simone-silvestri/TimestepperTestCases.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/simone-silvestri/TimestepperTestCases.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/simone-silvestri/TimestepperTestCases.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/simone-silvestri/TimestepperTestCases.jl)

**TimestepperTestCases.jl** is a Julia package that provides test cases and diagnostics for evaluating time-stepping schemes in free-surface ocean models. This repository complements the paper **"A Low-Storage Runge-Kutta Framework for Non-Linear Free-Surface Ocean Models"** by providing:

- Three idealized test cases spanning linear to fully turbulent regimes
- Exact variance budget diagnostics for quantifying numerical mixing
- Linear stability analysis tools
- Notebooks for reproducing figures and analysis from the paper

## Overview

This package implements test cases designed to highlight differences in performance between time discretization schemes, particularly comparing low-storage Runge-Kutta (RK) methods with Adams-Bashforth (AB2) schemes. The test cases demonstrate that RK3 time-stepping substantially reduces numerical mixing compared to AB2, preserving stratification and mesoscale variability in climate-relevant configurations.

### Key Findings

- **RK3 reduces numerical mixing**: Physical-to-numerical mixing ratios exceed 2 throughout most of the water column for RK3, compared to near-unity for AB2
- **Temporal discretization matters**: In linear regimes, temporal discretization effects can rival those from low-order spatial schemes
- **Stability constraints**: Split-explicit RK implementations require the averaged barotropic transport velocity for tracer advection to maintain baroclinic stability

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/simone-silvestri/TimestepperTestCases.jl")
```

## Quick Start

```julia
using TimestepperTestCases

# Run internal tide test case with RK3
sim = internal_tide(:SplitRungeKutta3)

# Run idealized coast test case
sim = idealized_coast(:SplitRungeKutta3; forced=true)

# Run channel simulation
sim = channel_simulation(timestepper=:SplitRungeKutta3)
```

See the [documentation](https://simone-silvestri.github.io/TimestepperTestCases.jl/stable/) for detailed usage instructions.

## Repository Structure

```
TimestepperTestCases.jl/
├── src/                                    # Source code
│   ├── TimestepperTestCases.jl           # Main module
│   ├── internal_tide.jl                  # Internal tide test case
│   ├── idealized_coast.jl                # Idealized coast test case
│   ├── channel_simulation.jl             # Channel simulation
│   ├── diagnostics.jl                    # RPE/APE diagnostics
│   ├── load_internal_tide_case.jl        # Load internal tide output
│   ├── load_idealized_coast_case.jl      # Load idealized coast output
│   ├── load_channel_case.jl              # Load channel output
│   └── BuoyancyVarianceDissipationComputations/  # Buoyancy variance diagnostics
│       ├── BuoyancyVarianceDissipationComputations.jl
│       ├── compute_dissipation.jl
│       ├── update_fluxes.jl
│       ├── advective_dissipation.jl
│       └── flatten_dissipation_fields.jl
├── notebooks/                             # Jupyter notebooks
│   ├── figure_1.ipynb                    # Amplification factors
│   ├── figure_2.ipynb                    # Linear stability
│   ├── linear_stability.ipynb            # Stability analysis
│   ├── Results_internal_tide.ipynb      # Internal tide results
│   ├── Results_idealized_coast.ipynb    # Idealized coast results
│   └── Results_channel.ipynb            # Channel results
├── docs/                                  # Documentation
│   ├── src/                              # Documentation source
│   │   ├── index.md                     # Home page
│   │   ├── getting_started.md           # Installation and usage
│   │   ├── experiments/                 # Test case descriptions
│   │   │   ├── internal_tide.md
│   │   │   ├── idealized_coast.md
│   │   │   └── channel.md
│   │   ├── diagnostics.md               # Variance budget methodology
│   │   ├── stability.md                 # Stability analysis
│   │   ├── notebooks.md                 # Notebook descriptions
│   │   └── api_reference.md            # API documentation
│   └── make.jl                          # Documentation build script
├── test/                                  # Test suite
│   └── runtests.jl
└── README.md                             # This file
```

## Experiments

### 1. Internal Tide

A linear test case simulating tidal flow over a Gaussian seamount. This configuration isolates temporal discretization effects, as spatial advection plays a secondary role.

**Key features:**
- 2 km deep domain with 250 m tall seamount
- Oscillatory tidal forcing (M₂ tide)
- Initial stratification `N² = 10⁻⁴ s⁻²`
- 40-day simulation

**Results:** RK3 schemes retain more energy, preserve sharper structures, and show lower numerical diffusivity than AB2.

### 2. Idealized Coastal Baroclinic Adjustment

A nonlinear test case with sloping bottom and initial meridional salinity gradient that drives baroclinic instabilities.

**Key features:**
- 192 km × 192 km × 103 m domain
- Linearly sloping bathymetry
- Meridional salinity gradient
- Optional wind forcing
- 40-day simulation

**Results:** RK3-SE maintains sharper stratification and stronger submesoscale variability, while AB2-SE suppresses instabilities through excessive numerical mixing.

### 3. Equilibrated Channel Flow

A long-term equilibrated configuration showing interaction between numerical and physical mixing.

**Key features:**
- 1000 km × 2000 km × 3 km periodic channel
- Sinusoidal wind stress and variable heat flux
- CATKE closure with background diffusivity
- 40-year simulation (last 5 years averaged)

**Results:** RK3-SE achieves `κ_phys / κ_num > 2` throughout most of the water column, while AB2-SE hovers around unity, indicating numerical mixing can overwhelm physical processes.

## Diagnostics

The package provides comprehensive diagnostics for quantifying numerical mixing:

- **Variance Budgets**: Exact, locally computed dissipation rates
- **Numerical Diffusivity**: Effective diffusivity from discretization errors
- **Reference Potential Energy (RPE)**: Quantifies irreversible mixing
- **Available Potential Energy (APE)**: Energy available for instabilities
- **Kinetic Energy**: Total and mean kinetic energy time series

See the [diagnostics documentation](https://simone-silvestri.github.io/TimestepperTestCases.jl/stable/diagnostics/) for details.

## Notebooks

Jupyter notebooks are provided for:
- Reproducing figures from the paper
- Analyzing simulation results
- Performing linear stability analysis
- Computing custom diagnostics

All notebooks are located in the `notebooks/` directory. See the [notebooks documentation](https://simone-silvestri.github.io/TimestepperTestCases.jl/stable/notebooks/) for descriptions and usage.

## Documentation

Full documentation is available at:
- **Stable**: https://simone-silvestri.github.io/TimestepperTestCases.jl/stable/
- **Dev**: https://simone-silvestri.github.io/TimestepperTestCases.jl/dev/

Build documentation locally:

```bash
julia --project=docs docs/make.jl
```

## Citation

If you use this package in your research, please cite:

```bibtex
@article{silvestri2025timestepper,
  title={A Low-Storage Runge-Kutta Framework for Non-Linear Free-Surface Ocean Models},
  author={Silvestri, Simone and others},
  journal={Journal of Advances in Modeling Earth Systems},
  year={2025}
}
```

## Related Packages

- [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) - The ocean model framework used by this package

## License

This package is licensed under the MIT License. See `LICENSE` for details.
