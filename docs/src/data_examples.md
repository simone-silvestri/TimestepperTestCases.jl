# Reproducing the paper

The simulation outputs behind the manuscript are archived in Zenodo record 24315107
([doi:10.5281/zenodo.24315107](https://doi.org/10.5281/zenodo.24315107), CC BY 4.0). The record holds three campaigns as
individual JLD2 files, `internal_tide` (12 GB), `channel` (25 GB) and `global_ocean` (100 GB), together with SHA-256
checksums, a file manifest, the final global timings, and a snapshot of the analysis source. The theory figures need no
data: they are computed directly by `TimestepperTheory`.

## Figure map

| Manuscript | Script | Data | Rendered page |
| --- | --- | --- | --- |
| Figures 1–9 | `lib/TimestepperTheory/examples/time_discretization.jl` | none | [Manuscript example](generated/time_discretization.md) |
| Figures 10–13 | `examples/figures/internal_tide_figures.jl` | `internal_tide` | [Internal tide figures](generated/internal_tide_figures.md) |
| Figures 14–15 | `examples/figures/channel.jl` | `channel` (averages files) | [Channel figures](generated/channel.md) |
| Figure 16 | `examples/figures/global_ocean_results.jl` | `global_ocean` (dissipation files) | [Global numerical diffusivity](generated/global_ocean_results.md) |
| Figure 17 | `examples/figures/global_surface_results.jl` | `global_ocean` (surface files) | [Global surface velocity](generated/global_surface_results.md) |
| Table 3 | `examples/figures/global_ocean.jl` | `global_ocean/reported_timings.csv` | [Global computational cost](generated/global_ocean.md) |
| Figure 19 | `lib/TimestepperTheory/examples/rotation.jl` | none | — |

The global surface tracer histories ([Surface tracers](generated/global_surface_averages.md)) use the `global_ocean` average
files and do not appear in the manuscript. The figures of the rendered pages are the published ones, stored in
`docs/figures/`; the global figures were redrawn from the downloaded record and match the published versions pixel for
pixel.

## Download

From the repository root, download and verify the campaigns you need:

```sh
julia scripts/download_data.jl /path/to/data internal_tide channel global_ocean
```

Without campaign names, the script downloads all three; a numeric first argument selects another version of the record.
It uses Julia standard libraries only, restores the `internal_tide/`, `channel/` and `global_ocean/` folders, and verifies
every file against `SHA256SUMS`. Files already present with a matching checksum are reused, and an interrupted transfer
resumes from its `.partial` file, both within a run and when the script is rerun. Version 1.0 of `SHA256SUMS` also lists
`global_ocean/global_ocean_SRK3-SE_cost.jld2`, which the record does not provide; the script skips it, and no figure uses
it. The checksums and download locations are also listed in `examples/datasets.toml`.

## Run the figure scripts

Instantiate the examples environment once, then pass each campaign directory to its script:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/figures/internal_tide_figures.jl /path/to/data/internal_tide
julia --project=examples examples/figures/channel.jl /path/to/data/channel
julia --project=examples examples/figures/global_ocean_results.jl /path/to/data/global_ocean
julia --project=examples examples/figures/global_surface_results.jl /path/to/data/global_ocean
julia --project=examples examples/figures/global_surface_averages.jl /path/to/data/global_ocean
julia --project=examples examples/figures/global_ocean.jl
```

Without a directory argument, `internal_tide_figures.jl` and the global scripts read `internal_tide/` and `global_ocean/` in
the repository root, the layout produced by `julia scripts/download_data.jl . internal_tide global_ocean`. `channel.jl`
and `global_ocean.jl` download a missing file from the record into `data/<campaign>/` and verify its checksum. The
internal-tide figures are written to `<input directory>/figures/`; the channel and global figures, as PNG and PDF, to
`examples/output/`. Both locations are ignored by Git, and the original data files are never modified.

## Reductions

The scripts reduce the three-dimensional outputs once and cache the result, so that later runs redraw the figures in
seconds. Halos are removed and face values are interpolated to cell centers before any reduction.

- **Channel:** final-window zonal means of `Abx`, `Aby`, `Abz`, `Gbx`, `Gby`, `Gbz` and `Dbz`, with `y`, `z` and the output
  time, saved as `data/channel_simulation/compact_averages_<label>.jld2`.
- **Global diffusivity:** the variance tendency P and the squared gradient G of every complete 90-day window, summed by
  depth, by distance from topography and by latitude band, saved as `data/global_ocean/compact_global_ocean_<label>_dissipation.jld2`.
  The diffusivity κ = −ΣP/(2ΣG) is formed only after summation. The six-case comparison uses the eight complete windows
  ending at days 1170–1800.
- **Global surface:** the final daily snapshot of surface kinetic energy and the regional velocity spectra of 74 snapshots
  at five-day spacing over the final year, saved next to the input as `global_ocean_<label>_velocity_diagnostics.jld2`.
- **Global surface tracers:** area-weighted surface temperature and salinity of each output window, saved as
  `data/global_ocean/compact_global_ocean_<label>_average.jld2`.

To redraw from a cache, pass the directory that holds it; delete a cache file to recompute it.

## Contents and limits of the record

- WRK3-SE-SM05 in file names is labelled WRK3-SE-F1 in the figures.
- The global `_surface` files are a subset of the original output: daily `u` and `v` over the final 365 days, with the
  original grid, iterations, times and field metadata, compressed with shuffle and Zstandard. The `_average` and
  `_dissipation` files are the original outputs.
- `reported_timings.csv` holds the final timings of Table 3. The raw timing logs and the per-run cost records are not in
  the record.
- The source snapshot of version 1.0 (`analysis_source.tar.gz`) sets the channel sponge width to 9Ly/10 and leaves the
  global run settings (time steps, duration, budget windows, cold start, drag and vertical Courant limits) to the
  campaign driver. The channel runs used Ly/10, and the code in this repository sets the global run settings as the
  defaults of `global_ocean`.

## Documentation figures

The documentation build does not read simulation data. For each data-driven page it copies the figure from
`examples/output/`, when a local run has produced one, or otherwise from the committed copy in `docs/figures/`. On CI a
missing figure stops the build. After regenerating a figure, copy the PNG to `docs/figures/` so that it reaches the site.

`scripts/archive_global_surface.jl` prepares the global surface subset for a new version of the record: it copies the
selected arrays without changing precision, verifies every copied array and checks that each file loads as a
`FieldTimeSeries`.
