# Data-driven Julia examples

The examples read local campaign outputs or compact JLD2 datasets. Download locations in `examples/datasets.toml` are placeholders until the datasets are published. A missing URL produces an error asking for a local directory; it does not substitute synthetic results.

## Local execution

From the repository root:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/internal_tide.jl /path/to/internal_tide
julia --project=examples examples/channel.jl       /path/to/channel_simulation
julia --project=examples examples/global_ocean.jl  /path/to/global_ocean
```

Omitting the directory uses the corresponding case directory in this checkout. Local data always take precedence over downloads. The global example accepts a list of case labels after the directory for partial timing datasets; include AB2-SE because it defines the relative throughput.

```sh
julia --project=examples examples/global_ocean.jl /path/to/global_ocean AB2-SE WRK3-SE-SM05
```

## Computations and outputs

- Internal tide: numerical diffusivity from the ratio of horizontally averaged variance destruction and squared gradients, depth and temporal averages, and depth–time differences against WRK3-SE.
- Channel: numerical diffusivity profiles, physical-to-numerical diffusivity ratios, and signed variance-tendency differences against WRK3-SE. The final saved averaging window is used.
- Global: simulated years per wall-clock day and throughput relative to AB2-SE, calculated using the nominal production timestep.

Computations are shown in the rendered examples. Plotting helpers are in `examples/plotting.jl`. Figures are written to `examples/output/` as PNG and PDF; this directory is ignored by Git. These examples cover the diagnostics listed above; internal-tide snapshots, rectified-flow convergence, energetics, and global surface diagnostics remain in the locally retained notebooks pending migration and validation.

## Compact datasets

When reading the original internal-tide or channel outputs, `examples/compact_data.jl` extracts the required arrays without constructing an Oceananigans grid. This supports the older serialized grid layouts in the campaign files. Halos are removed and face values are interpolated to cell centers before reduction. No original data file is modified.

Compact files are saved in `data/internal_tide/` and `data/channel_simulation/`, with `compact_` prepended to the original filename. To use a compact dataset directly, pass that directory to the example. Channel files contain the final-window zonal means of Abx, Aby, Abz, Gbx, Gby, Gbz, and Dbz, plus y, z, and the output time. Internal-tide files contain centered Abx, Abz, and Gbz on x–z–time coordinates, the wet mask, and vertical cell widths. Global cost records are already small and are used directly.

## Publication

Upload the compact files, a description of the selected averaging windows and units, the simulation code revision, and the case parameters. Fill in each dataset's `doi`, version-specific `base_url`, and per-file SHA-256 hash in `examples/datasets.toml`. The hash table uses quoted filenames as keys. Downloaded files are checked against those hashes before being placed in the ignored `data/` cache.

```toml
[channel_simulation]
doi = "" # Add the published DOI.
base_url = "" # Add the version-specific files URL.
[channel_simulation.sha256]
# "compact_averages_AB2-SE.jld2" = "<SHA-256 of the uploaded file>"
```

The documentation build renders the example source and includes figures that have been generated locally. It does not download unpublished datasets or run large simulations.
