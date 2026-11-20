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
│   ├── figures/                 # Data analysis and figure drivers
│   ├── shared/                  # Data loading and shared plotting helpers
│   └── datasets.toml            # Dataset download locations and checksums
├── scripts/                     # Archive preparation and verified downloads
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
julia --project=examples examples/figures/internal_tide.jl /path/to/internal_tide
julia --project=examples examples/figures/channel.jl /path/to/channel_simulation
julia --project=examples examples/figures/global_ocean.jl /path/to/global_ocean
julia --project=examples examples/figures/global_ocean_results.jl /path/to/global_ocean
```

See the [data guide](docs/src/data_examples.md) for commands and diagnostic coverage. Locally retained notebooks are ignored and excluded from the published tree.

## Reproduce the manuscript figures

The simulation output is archived on Zenodo (Silvestri et al., 2026, [doi:10.5281/zenodo.24315107](https://doi.org/10.5281/zenodo.24315107), CC BY 4.0). `scripts/download_data.jl` downloads the selected campaigns of record 24315107, restores their folders, and verifies every file against the record's SHA-256 checksums with Julia standard libraries only. The record holds about 137 GB, so download only the campaigns you need:

| Campaign | Download | Figure driver | Manuscript figures |
|---|---|---|---|
| — | — | `lib/TimestepperTheory/examples/time_discretization.jl` | 1–9 |
| — | — | `lib/TimestepperTheory/examples/rotation.jl` | 19 |
| `internal_tide` | 12 GB | `examples/figures/internal_tide_figures.jl` | 10–13 |
| `channel` | 25 GB | `examples/figures/channel.jl` | 14–15 |
| `global_ocean` | 100 GB | `examples/figures/global_ocean_results.jl`, `examples/figures/global_surface_results.jl` | 16–17 |

From the repository root, after instantiating the examples environment as above:

```sh
julia scripts/download_data.jl /path/to/data internal_tide channel global_ocean
julia --project=examples examples/figures/internal_tide_figures.jl /path/to/data/internal_tide
julia --project=examples examples/figures/channel.jl /path/to/data/channel
julia --project=examples examples/figures/global_ocean_results.jl /path/to/data/global_ocean
julia --project=examples examples/figures/global_surface_results.jl /path/to/data/global_ocean
julia --project=lib/TimestepperTheory/examples --threads=auto lib/TimestepperTheory/examples/rotation.jl
```

The internal-tide figures are written to `/path/to/data/internal_tide/figures/`, the channel and global figures to `examples/output/`, and the theory figures to `lib/TimestepperTheory/output/`. The global timings of the manuscript's cost table are in `global_ocean/reported_timings.csv`; the record does not contain the raw timing logs.

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

The documentation workflow builds the site for pull requests and publishes successful `main` builds to GitHub Pages. In the repository settings, select **Pages → Build and deployment → Source → GitHub Actions**. Each build also uploads a downloadable `documentation` artifact.

## Tests

```sh
julia --project=lib/TimestepperTheory -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

The theory tests check modal identities, filter moments, stage updates, amplification factors, and temporal convergence. The root tests check baroclinic stability limits and require the Oceananigans environment.

## Citation

If you use this repository, please cite the accompanying manuscript, and the data record when you use the archived simulation output:

```bibtex
@unpublished{silvestri2026conservative,
  title  = {A conservative {Runge--Kutta} framework for nonlinear free-surface ocean models},
  author = {Silvestri, Simone and Campin, Jean-Michel and Wagner, Gregory and Constantinou, Navid C. and Bishnu, Siddhartha and Lee, Xin Kai and Ferrari, Raffaele},
  year   = {2026},
  note   = {Manuscript}
}

@dataset{silvestri2026data,
  title     = {Reproducibility data for ``A conservative {Runge--Kutta} framework for nonlinear free-surface ocean models''},
  author    = {Silvestri, Simone},
  year      = {2026},
  publisher = {Zenodo},
  version   = {1.0},
  doi       = {10.5281/zenodo.24315107}
}
```

## License

See [LICENSE](LICENSE).
