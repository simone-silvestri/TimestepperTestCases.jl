# Data-driven Julia examples

The examples read local campaign outputs or compact JLD2 datasets. Download locations in `examples/datasets.toml` are placeholders until the datasets are published. A missing URL produces an error asking for a local directory; it does not substitute synthetic results.

## Local execution

From the repository root:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/internal_tide_figures.jl /path/to/internal_tide
julia --project=examples examples/channel.jl       /path/to/channel_simulation
julia --project=examples examples/global_ocean.jl  /path/to/global_ocean
julia --project=examples examples/global_ocean_results.jl /path/to/global_ocean/raw
julia --project=examples examples/global_surface_averages.jl /path/to/global_ocean/global
```

Omitting the directory uses the corresponding case directory in this checkout: `internal_tide/` for the manuscript internal-tide campaign, `channel_simulation/` for the channel, and `global_ocean/` for the global timings and results. Local data always take precedence over downloads. The global example accepts a list of case labels after the directory for partial timing datasets; include AB2-SE because it defines the relative throughput.

```sh
julia --project=examples examples/global_ocean.jl /path/to/global_ocean AB2-SE WRK3-SE-SM05
```

## Computations and outputs

- Internal tide: numerical diffusivity from the ratio of horizontally averaged variance destruction and squared gradients, depth and temporal averages, and depth–time differences against WRK3-SE.
- Channel: numerical diffusivity profiles, physical-to-numerical diffusivity ratios, and signed variance-tendency differences against WRK3-SE. The final saved averaging window is used.
- Global cost: simulated years per wall-clock day and throughput relative to AB2-SE, calculated using the nominal production timestep.
- Global results: numerical diffusivity profiles, profiles grouped by distance from topography, latitude–depth sections, and whole-ocean and deep-ocean histories from the buoyancy-variance budgets.

Computations are shown in the rendered examples. Plotting helpers are in `examples/plotting.jl`. Global and channel figures are written to `examples/output/` as PNG and PDF; this directory is ignored by Git. The complete internal-tide renderer also includes snapshots, rectified-flow convergence, and energetics.

## Compact datasets

When reading original outputs, `examples/compact_data.jl` reads the stored grid using the pinned Oceananigans environment and extracts the required arrays. Halos are removed and face values are interpolated to cell centers before reduction. No original data file is modified.

Compact files are saved in `data/internal_tide_weno_new/` and `data/channel_simulation/`, with `compact_` prepended to the original filename. To use a compact dataset directly, pass that directory to the example. Channel files contain the final-window zonal means of Abx, Aby, Abz, Gbx, Gby, Gbz, and Dbz, plus y, z, and the output time. Internal-tide files contain centered Abx, Abz, and Gbz on x–z–time coordinates, the wet mask, and vertical cell widths. Global cost records are already small and are used directly.

The global diffusivity example reads `global_ocean_<label>_dissipation.jld2` and saves compact reductions in `data/global_ocean/`. Each reduction retains the variance tendency and squared gradient separately for every complete averaging window, depth, distance from topography, and latitude band. Pass the compact directory to redraw without reading the original files. The six-case comparison uses the eight complete 90-day windows ending at days 1170–1800.

`examples/internal_tide_figures.jl` renders all five internal-tide figures, including the snapshots, energetics, diffusivity, and timestep-refinement comparison. It reads the thirteen completed runs and writes PDF and PNG files to `<input directory>/figures/`. The original output files are unchanged.

The global campaign ran with `statistics=false`. Its `_average` files contain surface averages; volume-integrated kinetic energy, temperature, and salinity histories were not saved and cannot be recovered from these outputs. The diffusivity figures use the saved buoyancy-variance budgets. `global_surface_averages.jl` computes area-weighted surface temperature and salinity from the `_average` files and compares them with AB2-SE at matching output times.

## Publication

Upload the compact files, a description of the selected averaging windows and units, the simulation code revision, and the case parameters. Fill in each dataset's `doi`, version-specific `base_url`, and per-file SHA-256 hash in `examples/datasets.toml`. The hash table uses quoted filenames as keys. Downloaded files are checked against those hashes before being placed in the ignored `data/` cache.

```toml
[channel_simulation]
doi = "" # Add the published DOI.
base_url = "" # Add the version-specific files URL.
[channel_simulation.sha256]
# "compact_averages_AB2-SE.jld2" = "<SHA-256 of the uploaded file>"
```

The documentation build renders the example source and includes figures that have been generated locally. Where `examples/output/` holds no such figure, it uses the copy in `docs/figures/`, which is how the global results reach the deployed site; copy a figure there again after regenerating it. The build does not download unpublished datasets or run large simulations.

For campaign timings, use the same discretization as the simulation. For example, profile
`global_ocean_cost(discretization("WRK3-UP"); arch=GPU())` to measure third-order upwind tracer advection;
`WRK3-SE` uses WENO9. For an implicit-solver sensitivity run, pass
`free_surface=ImplicitFreeSurface(solver_method=:PreconditionedConjugateGradient; reltol=1e-8, maxiter=1000)`
to either `global_ocean` or `global_ocean_cost` for `WRK3-IM`. The progress log reports the iteration count and residual.
