# Getting started

## Installation

The package is not registered. Clone the repository and work from its root:

```sh
git clone https://github.com/simone-silvestri/TimestepperTestCases.jl
cd TimestepperTestCases.jl
```

The repository holds four Julia environments, each pinned by its own manifest. Use Julia 1.12 or newer.

| Environment | Purpose | Set up with |
| --- | --- | --- |
| `.` | Simulations and diagnostics (`TimestepperTestCases`, pinned Oceananigans revision) | `julia --project=. -e 'using Pkg; Pkg.instantiate()'` |
| `examples` | Figure scripts that read the published data | `julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'` |
| `lib/TimestepperTheory/examples` | Time-discretization theory and its figures | `julia --project=lib/TimestepperTheory/examples -e 'using Pkg; Pkg.develop(path="lib/TimestepperTheory"); Pkg.instantiate()'` |
| `docs` | This documentation | `julia --project=docs -e 'using Pkg; Pkg.develop([PackageSpec(path=pwd()), PackageSpec(path="lib/TimestepperTheory")]); Pkg.instantiate()'` |

## Discretizations

Every experiment runs the same set of time discretizations, identified by a label that combines the outer scheme, the
barotropic solver and, where relevant, the averaging filter or the tracer advection:

```julia
using TimestepperTestCases

[d.label for d in discretizations()]   # AB2-SE, WRK3-SE, WRK3-SE-SM05, SRK3-SE, MRK4-SE, WRK3-IM, WRK3-UP
d = discretization("WRK3-SE")
```

The figures label WRK3-SE-SM05 as WRK3-SE-F1. Each experiment derives the timestep and the free surface of a label from
its own stability parameters.

## Run an experiment

```julia
using TimestepperTestCases

internal_tide(discretization("WRK3-SE"))          # 40-day internal tide on the CPU
channel_simulation(discretization("WRK3-SE"))     # re-entrant channel, GPU by default
```

The quarter-degree global configuration lives in a package extension that loads with NumericalEarth and
WorldOceanAtlasTools and runs on a GPU. The campaign drivers in `examples/simulations/` run every discretization of an
experiment in turn. The experiment pages describe each configuration, its outputs and its manuscript figures:

- [Internal tide](experiments/internal_tide.md)
- [Re-entrant channel](experiments/channel.md)
- [Quarter-degree global ocean](experiments/global_ocean.md)

## Reproduce the figures

The outputs of all three campaigns are published on Zenodo. [Reproducing the paper](data_examples.md) explains how to
download them and which script draws each manuscript figure. The [variance-budget diagnostics](diagnostics.md) page
defines the numerical diffusivity used throughout, and the [time-discretization theory](theory.md) pages cover Sections
2–4 of the manuscript.
