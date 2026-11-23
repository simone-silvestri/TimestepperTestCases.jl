# Quarter-degree global ocean

The global configuration tests the time discretizations in a realistic, eddy-permitting ocean. We run the eORCA025
tripolar grid for five years, starting from ECCO2 and forcing the ocean only with an idealized zonal wind stress and a
restoring of surface temperature and salinity. Every discretization runs with the same physics, so the runs differ only
in how they step in time. We then compare their numerical diffusivity, diagnosed from the buoyancy-variance budget
([diagnostics](../diagnostics.md)), their surface velocity, and their cost per simulated year.

## What you need

The global case lives in a package extension, `TimestepperTestCasesNumericalEarthExt`, which Julia loads once
NumericalEarth and WorldOceanAtlasTools are loaded alongside TimestepperTestCases. Both are optional dependencies of the
package and have to be installed in the active environment. The runs are meant for a GPU; the published timings were
measured on one Nvidia H100. On first use NumericalEarth downloads the eORCA025 mesh and bathymetry, the ECCO2 state of
1 January 1993, and the monthly World Ocean Atlas climatologies used for restoring.

## The configuration

| | |
| :-- | :-- |
| Grid | eORCA025 tripolar mesh, largest connected basin only, columns shallower than 30 m treated as land |
| Vertical grid | 100 levels over 5000 m, stretched to a 5 m top cell, z-star coordinate |
| Initial state | ECCO2 temperature and salinity of 1 January 1993, ocean at rest |
| Wind stress | steady and zonal, piecewise cubic in latitude, between −0.1 and 0.2 N m⁻², zero poleward of 70° |
| Surface restoring | temperature and salinity towards the monthly World Ocean Atlas, piston velocity 1/6 m day⁻¹ (salinity through a virtual salt flux) |
| Equation of state | linear, thermal expansion 1.67 × 10⁻⁴ K⁻¹, haline contraction 7.80 × 10⁻⁴ psu⁻¹ |
| Vertical mixing | CATKE; no explicit horizontal viscosity or diffusivity |
| Bottom drag | quadratic, coefficient 0.003 |
| Momentum advection | vector-invariant WENO, ninth order for vorticity and fifth order for the other terms |
| Tracer advection | ninth-order WENO; third-order upwind for WRK3-UP |
| Vertical advection | adaptively implicit above a vertical Courant number of 1 (0.5 for AB2-SE), for momentum and tracers |

There is no atmosphere, no freshwater flux and no sea ice, so the restoring is the only surface buoyancy forcing.

## Time steps and spin-up

The three-stage Runge–Kutta schemes step at 45 minutes and MRK4-SE at 72 minutes. AB2-SE runs at 11 minutes: the step
implied by its stability ratio, 13 minutes, is unstable in this configuration.

| Case | Outer scheme | Barotropic solver | Δt |
| :-- | :-- | :-- | --: |
| AB2-SE | quasi-AB2 | forward–backward substeps, SM05 filter | 11 min |
| WRK3-SE-SM05 | WRK3 | forward–backward substeps, SM05 filter | 45 min |
| WRK3-SE | WRK3 | RK3 substeps, optimized filter | 45 min |
| SRK3-SE | SRK3 | RK3 substeps, optimized filter | 45 min |
| MRK4-SE | MRK4 | RK3 substeps, optimized filter | 72 min |
| WRK3-UP | WRK3 | RK3 substeps, optimized filter | 45 min |
| WRK3-IM | WRK3 | implicit | 45 min |

WRK3-SE-SM05 is labelled WRK3-SE-F1 in the manuscript. WRK3-IM is only timed; the other six cases run for five years.
The number of barotropic substeps is fixed for each case: it is the smallest multiple of 8 that keeps the substep Courant
number below 0.8 for forward–backward and 1.3 for RK3 substeps, measured on the fastest gravity wave of the grid.

Starting from rest, every case first runs for 360 days at a third of its time step. AB2-SE needs a gentler start: it runs
at 1.741 minutes for the first 360 days and at 5.5 minutes up to day 960 before switching to 11 minutes. The variance
budget is switched on at day 1095 and averaged over 90-day windows until the end of the run at day 1830.

## Running a case

With a GPU and the extension loaded, one case is a single call:

```julia
using TimestepperTestCases
using NumericalEarth, WorldOceanAtlasTools
using Oceananigans
using CUDA

arch = GPU()
grid = global_ocean_grid(arch)

result = global_ocean(discretization("WRK3-SE"); arch, grid)
```

The defaults reproduce the published runs. `result.ocean` is the Oceananigans `Simulation`, and `result` also carries the
wall time, the number of iterations, the time per step and the simulated years per day. Every setting above is a keyword
of `global_ocean`; the ones that shape the run are

- `Δt`, `cold_start_Δt`, `cold_start_duration`, `handover_Δt` and `handover_duration` for the time steps of the spin-up;
- `stop_time`, `dissipation_start_time` and `dissipation_output_interval` for the length of the run and of the budget;
- `bottom_drag_coefficient`, `wind_stress_values` and `restoring_piston_velocity` for the forcing;
- `checkpoint_interval` and `pickup`, to write checkpoints and resume from the newest one.

`global_ocean_cost(discretization("SRK3-SE"); arch = GPU())` measures the time per step without writing any output: it
discards 200 warm-up steps and times the next 1000. `examples/simulations/global_ocean.jl` runs the whole campaign,
timing every case first and then integrating the six five-year cases.

## Output

Each run writes to the working directory, under the name `global_ocean_<label>`:

| File | Contents |
| :-- | :-- |
| `_surface.jld2` | daily surface `T`, `S`, `e`, `u`, `v` and `w` |
| `_average.jld2` | the same surface fields, averaged over 90-day windows |
| `_dissipation.jld2` | the three-dimensional variance tendencies `Abx`, `Aby`, `Abz` and squared buoyancy gradients `Gbx`, `Gby`, `Gbz`, averaged over 90-day windows from day 1095 |
| `_checkpoint_iteration<N>.jld2` | the newest checkpoint |
| `_cost.jld2` | the time steps, substeps and timings of the run |

`load_global_ocean(folder, label)` opens these files as `FieldTimeSeries`.

The published record ([doi:10.5281/zenodo.24315107](https://doi.org/10.5281/zenodo.24315107)) holds the `_average` and
`_dissipation` files of the six five-year cases, a `_surface` subset with the daily `u` and `v` of the final year, and the
timings of the manuscript in `reported_timings.csv`. To open it with `load_global_ocean`, pass `surface = false`, since
the published surface files hold the velocities only. [Reproducing the paper](../data_examples.md) explains how to
download it.

## Results

Over days 1080–1800 the numerical diffusivity of the Runge–Kutta cases stays close to that of AB2-SE. Only WRK3-UP
differs, with a global diffusivity about twice as large ([numerical diffusivity](../generated/global_ocean_results.md)).
Surface speed and regional kinetic-energy spectra are similar across cases, with less energy for WRK3-UP
([surface velocity](../generated/global_surface_results.md)), and the surface temperature and salinity evolve alike
([surface tracers](../generated/global_surface_averages.md)). The longer steps of the Runge–Kutta schemes more than make
up for their extra stages, while the implicit free surface halves the throughput
([computational cost](../generated/global_ocean.md)).
