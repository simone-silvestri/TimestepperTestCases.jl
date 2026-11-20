# Quarter-degree global ocean: configuration and run

The global configuration tests the time discretizations of the package in a realistic geometry: an eddy-permitting ocean
on the eORCA025 tripolar grid, initialized from ECCO2 and forced by an idealized zonal wind stress and by surface
restoring of temperature and salinity. Every discretization runs with the same physics, and the runs are compared through
the numerical diffusivity diagnosed from the buoyancy-variance budget (see [diagnostics](../diagnostics.md)), the surface
velocity and the cost per step. This page describes the configuration as the code builds it and how to run it; the
analysis pages are listed at the end.

## Requirements

- **Extension.** `global_ocean`, `global_ocean_grid`, `global_ocean_cost`, `global_ocean_vertical_discretization` and
  `global_ocean_barotropic_rate` have methods only in the package extension `TimestepperTestCasesNumericalEarthExt`.
  Julia loads it when NumericalEarth and WorldOceanAtlasTools are loaded next to TimestepperTestCases. Both are weak
  dependencies in `Project.toml` (NumericalEarth 0.8 or 0.9, WorldOceanAtlasTools 0.6.5) and have to be installed in the
  active environment.
- **Hardware.** The default architecture is `CPU()`; `examples/simulations/global_ocean.jl` runs on `GPU()` and loads
  CUDA. The published timings were measured on one Nvidia H100.
- **Input data**, downloaded on first use:
  - the eORCA025 mesh and bathymetry, fetched by NumericalEarth for `ORCAQuarter()`;
  - ECCO2 daily temperature and salinity of 1 January 1993 (`ECCO2Daily()`), the initial condition. Files missing from
    the NumericalEarth data directory are downloaded from the `data-v1` release of NumericalEarthArtifacts;
  - the monthly World Ocean Atlas temperature and salinity climatologies (`WOAMonthly()`), the restoring targets,
    downloaded by the World Ocean Atlas extension of NumericalEarth that WorldOceanAtlasTools activates.

## Configuration

| Component | Setting |
| :-- | :-- |
| Horizontal grid | eORCA025 tripolar mesh (`ORCAQuarter()`), periodic in longitude and folded at the northern boundary. Only the largest connected basin is kept (`major_basins = 1`), and columns shallower than 30 m are land (`minimum_depth`). Zero spacings and areas under the fold and in the southern halo are replaced by the nearest positive value on the same row. |
| Vertical grid | 100 levels over 5000 m, exponentially stretched to a 5 m top cell (`global_ocean_vertical_discretization`). The coordinate is mutable, so the model runs in z-star. |
| Bathymetry | eORCA025 bathymetry as a `GridFittedBottom` immersed boundary with an active-cells map |
| Initial conditions | ECCO2 temperature and salinity of 1 January 1993 (`init_date`), temperature clipped below at −1.8 °C; velocities and free surface at rest |
| Wind stress | Zonal and steady, piecewise cubic in latitude through the knots (−70, −45, −15, 0, 15, 45, 70)° with values (0, 0.2, −0.1, −0.02, −0.1, 0.1, 0) N m⁻², zero poleward of ±70°. It is divided by a reference density of 1026 kg m⁻³ and rotated onto the grid directions. |
| Surface restoring | Temperature and salinity restored towards the monthly WOA climatology with a piston velocity of 1/6 m day⁻¹, a 30-day time scale for the 5 m top cell. The restoring is a top flux (`SurfaceFluxRestoring`), so salinity is restored through a virtual salt flux and the interior budget carries no restoring term. |
| Other forcing | None: no atmosphere, no freshwater flux, no sea ice |
| Equation of state | `LinearEquationOfState()`, thermal expansion 1.67 × 10⁻⁴ K⁻¹ and haline contraction 7.80 × 10⁻⁴ psu⁻¹; reference density 1026 kg m⁻³ and heat capacity 3992 J kg⁻¹ K⁻¹ for the NumericalEarth flux conversions |
| Closure | CATKE only, vertically implicit, with Cᵇ = 0.01 in the mixing length and Cᵂϵ = 1 in the TKE equation; no explicit horizontal viscosity or diffusivity |
| Bottom drag | Quadratic, coefficient 0.003, background velocity 0.1 m s⁻¹ for the unforced barotropic tide |
| Momentum advection | Vector-invariant WENO: ninth order for vorticity, fifth order for divergence, kinetic-energy gradient and vertical advection, reducing to first-order upwind near boundaries (`split_momentum_advection`) |
| Tracer advection | Ninth-order WENO with stretched vertical reconstruction (`weno_tracer_advection`); third-order upwind for WRK3-UP |
| Vertical advection | Adaptive implicit (`AdaptiveVerticallyImplicitDiscretization`): the part of the vertical velocity above a vertical Courant number `cfl` is advected implicitly with first-order upwinding. `cfl = 1`, and 0.5 for AB2-SE. It applies to momentum and to the tracers, including the third-order upwind scheme of WRK3-UP. |
| Free surface | Split-explicit with a fixed substep count, implicit for WRK3-IM |
| Duration | 360-day cold start at Δt/3, then Δt up to `stop_time = 1830days` (five years). AB2-SE cold-starts at 1.741 minutes and runs at 5.5 minutes up to day 960 before its 11-minute step. The variance budget joins at day 1095 |

## Discretizations

`global_ocean(d::Discretization)` runs one row of `discretizations()`. The three-stage schemes step at 45 minutes,
MRK4-SE at 72 minutes and AB2-SE at 11 minutes (Table 3 of the manuscript): the step implied by the AB2 ratio, 13
minutes, is unstable in this configuration, and 11 minutes is the largest stable one. The barotropic substep count is
the smallest multiple of 8 that keeps the substep Courant number c₀kΔτ at or below 0.8 for FB and 1.3 for RK3 substeps, with c₀k the maximum over the
columns of √(gH) 2√(Δx⁻² + Δy⁻²) (`global_ocean_barotropic_rate`). It is sized for `nominal_timestep`.

| Label | `timestepper` | Free surface | Substeps | Averaging kernel | Slow forcing | Tracer advection | Δt |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- | --: |
| AB2-SE | `:QuasiAdamsBashforth2` | split-explicit | forward–backward | SM05 | frozen | WENO9 | 11 min |
| WRK3-SE-SM05 | `:SplitRungeKutta3` | split-explicit | forward–backward | SM05 | frozen | WENO9 | 45 min |
| WRK3-SE | `:SplitRungeKutta3` | split-explicit | RK3 | optimized asymmetric | progressive | WENO9 | 45 min |
| SRK3-SE | `:SSPRungeKutta3` | split-explicit | RK3 | optimized asymmetric | frozen | WENO9 | 45 min |
| MRK4-SE | `:ModifiedRungeKutta4` | split-explicit | RK3 | optimized asymmetric | progressive | WENO9 | 72 min |
| WRK3-IM | `:SplitRungeKutta3` | implicit | none | none | none | WENO9 | 45 min |
| WRK3-UP | `:SplitRungeKutta3` | split-explicit | RK3 | optimized asymmetric | progressive | upwind, third order | 45 min |

SM05 is `averaging_shape_function` and the optimized asymmetric kernel is `OptimizedAsymmetricAveragingKernel()`.
WRK3-SE-SM05 is labelled WRK3-SE-F1 in the manuscript. The published campaign integrated every row except WRK3-IM for
five years and timed WRK3-IM only; `global_ocean_discretizations()` returns the six integrated rows.

## Running

### One discretization

```julia
using TimestepperTestCases
using NumericalEarth, WorldOceanAtlasTools
using Oceananigans
using Oceananigans.Units
using CUDA

arch = GPU()
grid = global_ocean_grid(arch)

result = global_ocean(discretization("WRK3-SE"); arch, grid)
```

`result` is the named tuple `(; ocean, label, wall_time, iterations, seconds_per_step, nominal_timestep,
simulated_years_per_day)`, with `ocean` the Oceananigans `Simulation`. The run has four phases:

1. a cold start of `cold_start_duration = 360days` at `cold_start_Δt`, followed for AB2-SE by `handover_duration = 600days`
   at `handover_Δt`;
2. a production segment at `Δt` up to `cost_stop_time = min(stop_time, 365days)`, whose wall time per step is
   `seconds_per_step`;
3. if the model time is still below `dissipation_start_time` (day 1095 by default), a run up to it;
4. the buoyancy-variance budget and the statistics writer are attached, the budget is evaluated every iteration, and
   the run continues to `stop_time`.

`global_ocean(d::Discretization; arch, grid, kw...)` sets `free_surface`, `averaging_kernel`, `barotropic_timestepper`,
`slow_forcing`, `label` and, for WRK3-UP only, `tracer_advection` from `d`, and forwards every other keyword to
`global_ocean(timestepper::Symbol; kw...)`. The keywords that set the run are:

| Keyword | Default | Effect |
| :-- | :-- | :-- |
| `Δt` | from the scheme, see above | production time step |
| `nominal_timestep` | `Δt` | step the substeps are sized for and the throughput is computed with |
| `cold_start_Δt`, `cold_start_duration` | `Δt / 3` (1.741 minutes for AB2), `360days` | cold start |
| `handover_Δt`, `handover_duration` | `nothing` (5.5 minutes for AB2), `600days` | intermediate step between the cold start and `Δt` |
| `stop_time` | `1830days` | end of the run |
| `cost_stop_time` | `min(stop_time, 365days)` | end of the timed segment |
| `dissipation`, `dissipation_start_time` | `true`, `1095days` | buoyancy-variance budget and the time it is attached |
| `output`, `statistics` | `true`, `false` | surface, average and statistics writers |
| `surface_output_interval` | `TimeInterval(1days)` | surface snapshots |
| `dissipation_output_interval` | `AveragedTimeInterval(90days)` | averaging window of the `_average` and `_dissipation` files |
| `statistics_output_interval` | `AveragedTimeInterval(7days)` | averaging window of the `_statistics` file |
| `checkpoint_interval`, `pickup` | `60days`, `false` | checkpoints; `pickup = true` resumes from the newest one |
| `prefix`, `label` | `"global_ocean"`, from `d` | output file names `<prefix>_<label>_<kind>.jld2` |
| `closure`, `equation_of_state` | CATKE as above, `LinearEquationOfState()` | physics |
| `bottom_drag_coefficient`, `bottom_drag_background_velocity` | `0.003`, `0.1` | quadratic drag |
| `wind_stress`, `wind_stress_latitudes`, `wind_stress_values` | `true` and the knots above | wind forcing |
| `surface_restoring`, `restoring_piston_velocity`, `restoring_dataset` | `true`, `1/6`, `WOAMonthly()` | restoring |
| `init_date` | `DateTime(1993, 1, 1)` | ECCO2 initial date |

A resumed run (`pickup = true`) skips the cold start and the timed segment, so its `seconds_per_step` is `NaN`. The
defaults reproduce the published runs: 1830 days with 90-day averaging windows of the budget.

### Cost measurement

```julia
cost = global_ocean_cost(discretization("SRK3-SE"); arch = GPU())
cost.seconds_per_step, cost.simulated_years_per_day
```

`global_ocean_cost(discretization = :SplitRungeKutta3; arch = CPU(), grid = global_ocean_grid(arch), warmup_steps = 200,
measured_steps = 1000, kw...)` sizes the substeps for the nominal time step of the scheme and steps at a third of it,
the cold-start step: the work per step does not depend on the step, and the run stays far from its stability limit. The
`warmup_steps` absorb compilation and are not timed. The implicit free surface is measured at the nominal step, since
its solver iterations grow with the step. The run writes no fields and no budget; its record goes to
`global_ocean_cost_<label>_cost.jld2`. Throughput is Δt/(365 t_step) at the nominal step
(`global_ocean_simulated_years_per_day`).

### Campaign driver

`examples/simulations/global_ocean.jl` runs on the GPU. It first times all seven `discretizations()` with
`global_ocean_cost` and prints Δt, seconds per step and simulated years per day for each, then integrates every case of
`global_ocean_discretizations()` with `global_ocean(d; arch)` and the default keywords. From the repository root, in an
environment that provides NumericalEarth, WorldOceanAtlasTools and CUDA:

```sh
julia --project=/path/to/environment examples/simulations/global_ocean.jl
```

## Output

`global_ocean` writes into the working directory, with `<prefix>_<label>` as the base name (`global_ocean_WRK3-SE`, for
example):

| File | Contents | Schedule |
| :-- | :-- | :-- |
| `_surface.jld2` | `T`, `S`, `e`, `u`, `v`, `w` at the top level, with halos, in Float32 | `surface_output_interval`, from the start |
| `_average.jld2` | the same surface fields, time-averaged | `dissipation_output_interval`, from the start |
| `_dissipation.jld2` | three-dimensional `Abx`, `Aby`, `Abz` (advective variance production along each direction) and `Gbx`, `Gby`, `Gbz` (squared buoyancy gradients), all integrated over each cell volume and time-averaged, in Float32 | `dissipation_output_interval`, from `dissipation_start_time` |
| `_statistics.jld2` | time averages of `u`, `v`, `T`, `S`, `e`, `η`, `b`, their squares, and the CATKE dissipation χᴷ = 2κᶜ(∂z b)² V, in Float64 | `statistics_output_interval`, from `dissipation_start_time` |
| `_checkpoint_iteration<N>.jld2` | model state, newest only | `checkpoint_interval` |
| `_cost.jld2` | `Δt`, `nominal_timestep`, `substeps`, wall time of each phase, `iterations`, `production_steps`, `seconds_per_step`, `simulated_years_per_day`, `cost_stop_time`, `stop_time` | end of the run |

`load_global_ocean(folder, label)` opens these files as on-disk `FieldTimeSeries`: the budget terms under `:Abx` to
`:Gbz` with `:times`, the surface fields under `:T`, `:S`, `:u`, `:v`, `:w` with `:surface_times`, their averages under
`:T̄`, `:S̄`, `:ū`, `:v̄`, `:w̄` with `:average_times`, and the timings under `:cost`; `statistics = true` adds the
`_statistics` fields. `load_global_ocean_cases(folder; labels)` opens every label whose files are present. Its default
`labels` come from `global_ocean_discretizations()`, so WRK3-UP has to be passed explicitly.

### Published record

The `global_ocean` campaign of Zenodo record 24315107 ([doi:10.5281/zenodo.24315107](https://doi.org/10.5281/zenodo.24315107))
holds five-year runs (1830 days) of AB2-SE, WRK3-SE-SM05, WRK3-SE, SRK3-SE, MRK4-SE and WRK3-UP. For each it has:

- `global_ocean_<label>_average.jld2`, the surface averages as written;
- `global_ocean_<label>_dissipation.jld2`, the budget terms averaged over 90-day windows; the comparisons use the eight
  complete windows ending at days 1170–1800;
- `global_ocean_<label>_surface.jld2`, a subset of the surface output with daily `u` and `v` over the final 365 days.

`reported_timings.csv` holds the final timings of the cost table, WRK3-IM included. The record has no `_statistics` or
`_cost` files. `load_global_ocean` reads `T`, `S`, `u`, `v` and `w` from the `_surface` file, so on the published record
it is called with `surface = false`, which opens the budget and the averages; the published surface velocities are read
by the figure scripts. Download and figure commands are in [Reproducing the paper](../data_examples.md).

## Results

Over the budget window from day 1080 to day 1800, the numerical diffusivity of the Runge–Kutta configurations stays
close to AB2-SE, except for WRK3-UP, whose global value is about twice as large ([numerical
diffusivity](../generated/global_ocean_results.md)). Surface speed and regional kinetic-energy spectra are similar across
cases, with less spectral energy for WRK3-UP ([surface speed and spectra](../generated/global_surface_results.md));
surface temperature and salinity histories are compared in [surface averages](../generated/global_surface_averages.md).
The longer steps of the Runge–Kutta schemes more than compensate for their extra stages, and WRK3-IM reaches about half
the throughput of AB2-SE ([computational cost](../generated/global_ocean.md)).
