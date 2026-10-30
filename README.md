# TimestepperTestCases.jl

[![Documentation](https://img.shields.io/badge/docs-online-blue.svg)](https://simone-silvestri.github.io/TimestepperTestCases.jl/)

Julia implementations of ocean-model test cases, buoyancy-variance diagnostics, and the time-discretization analysis accompanying the manuscript *A conservative Runge–Kutta framework for nonlinear free-surface ocean models*.

## Packages

`TimestepperTestCases` in `src/` provides Oceananigans simulation setups and numerical-mixing diagnostics. [TimestepperTheory](lib/TimestepperTheory/README.md), in `lib/TimestepperTheory/`, provides amplification factors, averaging filters, modal evolution matrices, and direct and indirect coupling corrections.

The theory package runs on a CPU without Oceananigans or simulation output. Its Literate example reproduces the computations and figures of Sections 2–4. The simulation package includes internal tides, re-entrant channel flow, and an optional NumericalEarth extension for the global configuration.

## Repository layout

```text
TimestepperTestCases.jl/
├── src/                         # Simulation setups and variance diagnostics
├── ext/                         # NumericalEarth global configuration
├── lib/TimestepperTheory/       # Theory package, tests, and Literate example
├── examples/
│   ├── simulations/             # Internal-tide, channel, and global drivers
│   ├── internal_tide.jl         # Internal-tide diagnostics
│   ├── channel.jl               # Channel diagnostics
│   ├── global_ocean.jl          # Global throughput comparison
│   └── datasets.toml            # Dataset download locations and checksums
├── docs/                        # Documenter sources and build script
└── test/                        # Simulation-package tests
```

Local simulation data, generated figures, notebooks, and documentation builds are ignored by Git.

## Run the theory example

Run these commands from the repository root with Julia 1.12:

```sh
julia --project=lib/TimestepperTheory/examples -e 'using Pkg; Pkg.develop(path="lib/TimestepperTheory"); Pkg.instantiate()'
julia --project=lib/TimestepperTheory/examples lib/TimestepperTheory/examples/time_discretization.jl
```

The example saves nine figures and their numerical values as PDF, PNG, and CSV files in `lib/TimestepperTheory/output`. The computations appear in the Literate source; plotting helpers live in `time_discretization_figures.jl`. Build the [documentation](#documentation) to read the rendered example at `docs/build/generated/time_discretization.html`.

See the [example guide](docs/src/theory_example.md) for the experiments, parameter choices, and output files. Adding `--quick` to the computation command uses coarse parameter grids and overwrites the same output files; rerun without it to restore the manuscript grids.

## Run the ocean simulations

Instantiate the root environment, which specifies the Oceananigans revision and other simulation dependencies:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The simulation drivers are in `examples/simulations/`; run them from the repository root. The global tripolar configuration uses the NumericalEarth and WorldOceanAtlasTools extension.

## Analyze local or published data

The Julia examples accept a local data directory and generate figures in the ignored `examples/output/` directory:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/internal_tide.jl /path/to/internal_tide
julia --project=examples examples/channel.jl /path/to/channel_simulation
julia --project=examples examples/global_ocean.jl /path/to/global_ocean
```

Without a directory argument, each example looks in the corresponding case directory in this checkout. Internal-tide and channel outputs are reduced to compact plain-array JLD2 files for publication. Dataset URLs and checksums in `examples/datasets.toml` are placeholders pending upload. See the [data guide](docs/src/data_examples.md) for the covered diagnostics, compact formats, and publication steps. Locally retained notebooks are ignored and excluded from the published tree.

## Documentation

The Documenter site contains the simulation documentation, a guide to the theory example, its rendered computations and figures, and API references for both packages.

```sh
julia --project=docs -e 'using Pkg; Pkg.develop([PackageSpec(path=pwd()), PackageSpec(path="lib/TimestepperTheory")]); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Open `docs/build/index.html`. The default build executes the full theory example. To rebuild with existing figures:

```sh
julia --project=docs docs/make.jl --reuse-figures
```

The documentation workflow builds the site and uploads a `documentation` artifact. It does not publish the site.

## Tests

```sh
julia --project=lib/TimestepperTheory -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

The theory tests check modal identities, filter moments, stage updates, amplification factors, and temporal convergence. The root tests check baroclinic stability limits and require the Oceananigans environment.

## Citation

If you use this repository, please cite the accompanying manuscript:

```bibtex
@unpublished{silvestri2026conservative,
  title  = {A conservative {Runge--Kutta} framework for nonlinear free-surface ocean models},
  author = {Silvestri, Simone and Campin, Jean-Michel and Wagner, Gregory and Constantinou, Navid C. and Bishnu, Siddhartha and Lee, Xin Kai and Ferrari, Raffaele},
  year   = {2026},
  note   = {Manuscript}
}
```

## License

See [LICENSE](LICENSE).
