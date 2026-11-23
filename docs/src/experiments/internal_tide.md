# Internal tide: configuration and run

The internal-tide case drives a semidiurnal tide over a Gaussian seamount in a uniformly stratified, two-dimensional
(x, z) domain, following the internal-tide example of the Oceananigans documentation. The tidal excursion (2 km) is
small compared with the seamount width (20 km), so the flow is mostly linear and the differences between runs come
mainly from the time discretization. No closure is attached: all tracer and momentum dissipation is numerical. Every
run integrates 40 days and records the buoyancy-variance budget described in [diagnostics](../diagnostics.md), from
which the numerical diffusivity is computed. The [manuscript figures](../generated/internal_tide_figures.md) page
draws the results.

## Configuration

All values are set in `src/internal_tide.jl` (`internal_tide_parameters`, `internal_tide_grid`, `internal_tide`).

| Quantity | Value |
|:--|:--|
| Model | `HydrostaticFreeSurfaceModel`, `BuoyancyTracer()`, tracers `(:b, :c)` |
| Domain | `x ∈ [-1000, 1000]` km, periodic; `z ∈ [-2, 0]` km, bounded; topology `(Periodic, Flat, Bounded)` |
| Grid | `RectilinearGrid`, `Nx × Nz = 256 × 128`, `Δx = 7812.5` m, `Δz = 15.625` m at rest, halo `(6, 6)`, z-star `MutableVerticalDiscretization` (default `mutable = true`) |
| Topography | `GridFittedBottom` at `-H + h₀ exp(-x² / 2width²)`, `h₀ = 250` m, `width = 20` km |
| Rotation | `FPlane(f = -1.03126e-4)` s⁻¹ |
| Stratification | `b = Nᵢ² z`, `Nᵢ² = 10⁻⁴` s⁻² |
| Tidal forcing | body force on `u`: `A₂ sin(ω₂ t)`, `T₂ = 12.421` h, `ω₂ = 2π / T₂`, `U₂ = ϵ ω₂ width ≈ 0.281` m s⁻¹ with `ϵ = 0.1`, `A₂ = U₂ (ω₂² - f²) / ω₂ ≈ 1.82 × 10⁻⁵` m s⁻² |
| Initial state | `u = U₂`, `b = Nᵢ² z`; the passive tracer `c` is not set and stays zero |
| Closure | none: no viscosity, diffusivity or bottom drag |
| Tracer advection | `WENO(order = 7)` (`UpwindBiased(order = 3)` for WRK3-UP) |
| Momentum advection | `FluxFormAdvection` with `WENO(order = 5)` in x and z |
| Duration | `stop_time = 40days` |
| Architecture | CPU: `internal_tide_grid` uses the default architecture and `compute_tracer_dissipation!` launches on `CPU()` |

## Discretizations

The campaign runs the seven [`Discretization`](@ref)s returned by [`discretizations`](@ref) (`src/discretizations.jl`).
FB is `ForwardBackwardScheme()`, RK3 is `RungeKutta3Scheme()`; frozen and progressive refer to `FrozenSlowForcing()`
and `ProgressiveSlowForcing()`.

| Label | `timestepper` | Free surface, substeps | Averaging kernel | Slow forcing | Tracer advection | `Δt` |
|:--|:--|:--|:--|:--|:--|:--|
| AB2-SE | `:QuasiAdamsBashforth2` | split-explicit, FB, 32 | `averaging_shape_function` | frozen | WENO, 7th order | 196.3 s (3.27 min) |
| WRK3-SE-SM05 | `:SplitRungeKutta3` | split-explicit, FB, 96 | `averaging_shape_function` | frozen | WENO, 7th order | 676.6 s (11.28 min) |
| WRK3-SE | `:SplitRungeKutta3` | split-explicit, RK3, 64 | `OptimizedAsymmetricAveragingKernel` | progressive | WENO, 7th order | 676.6 s (11.28 min) |
| SRK3-SE | `:SSPRungeKutta3` | split-explicit, RK3, 64 | `OptimizedAsymmetricAveragingKernel` | frozen | WENO, 7th order | 676.6 s (11.28 min) |
| MRK4-SE | `:ModifiedRungeKutta4` | split-explicit, RK3, 96 | `OptimizedAsymmetricAveragingKernel` | progressive | WENO, 7th order | 1088.7 s (18.14 min) |
| WRK3-IM | `:SplitRungeKutta3` | `ImplicitFreeSurface()` | | | WENO, 7th order | 676.6 s (11.28 min) |
| WRK3-UP | `:SplitRungeKutta3` | split-explicit, RK3, 64 | `OptimizedAsymmetricAveragingKernel` | progressive | upwind, 3rd order | 676.6 s (11.28 min) |

The manuscript calls `averaging_shape_function` F1 (so WRK3-SE-SM05 appears there as WRK3-SE-F1) and
`OptimizedAsymmetricAveragingKernel` F3.

Neither `Δt` nor the substep count is stored. `timestep_and_free_surface` derives both from the named tuple
`internal_tide_stability_parameters() = (N² = 1e-4, H = 2000, Δx = 7812.5, safety = 1)`. [`baroclinic_timestep`](@ref)
gives `Δt = θ★ / (c₁ k)` with `c₁ = N H / π` and `k = π / Δx`, so that each scheme runs at its own imaginary-axis limit
`θ★` on the grid-scale first baroclinic mode: `θ★ = 0.5025` for quasi-AB2 (`χ = 0.1`), `√3` for the three-stage
schemes, and `2.787` for MRK4 with `a = 0.23`, the default of the `MRK4_A` environment variable read by
`campaign_timestepper`, so that MRK4 runs at 1.609 times the three-stage step. [`barotropic_substeps`](@ref) returns the smallest multiple of 8 that keeps the substep Courant
number `c₀ k Δτ`, with `c₀ = √(gH)`, at or below 0.8 for FB and 1.3 for RK3 substeps. In the timestep-refinement runs
the substep count is recomputed from the scaled step: 8 and 8 for AB2-SE, 16 and 8 for WRK3-SE and SRK3-SE, 24 and 16
for MRK4-SE, at `Δt / 4` and `Δt / 8` respectively.

## Running

`internal_tide` has two methods:

```julia
internal_tide(timestepper::Symbol;
              grid = internal_tide_grid(),
              barotropic_timestepper = ForwardBackwardScheme(),
              averaging_kernel = OptimizedAsymmetricAveragingKernel(),
              free_surface = SplitExplicitFreeSurface(grid; averaging_kernel,
                                                      timestepper = barotropic_timestepper,
                                                      substeps = internal_tide_substeps(barotropic_timestepper,
                                                                                        timestepper;
                                                                                        averaging_kernel)),
              free_surface_name = default_free_surface_name(free_surface),
              tracer_advection = nothing,
              momentum_advection = nothing,
              stop_time = 40days,
              output_schedule = TimeInterval(1hours),
              run_simulation = true,
              Δt = internal_tide_timestep(Val(timestepper)))

internal_tide(d::Discretization; grid = internal_tide_grid(), timestep_factor = 1, label = d.label, kw...)
```

The `Discretization` method computes `Δt` and the free surface with `timestep_and_free_surface`, multiplies `Δt` by
`timestep_factor` (recomputing the substeps when the factor differs from 1), and calls the symbol method with
`free_surface_name = label` and the tracer advection of `d`. The remaining keywords (`stop_time`, `output_schedule`,
`run_simulation`, `momentum_advection`) pass through. The symbol method on its own uses FB substeps with
`OptimizedAsymmetricAveragingKernel`, which matches none of the seven campaign rows; use the `Discretization` method to
reproduce a campaign run. Both methods return the `Simulation`, already run unless `run_simulation = false`.

```julia
using TimestepperTestCases
using Oceananigans
using Oceananigans.Units

# Writes internal_tide_SplitRungeKutta3_WRK3-SE.jld2 in the working directory
simulation = internal_tide(discretization("WRK3-SE"); stop_time = 2days)

# Build without running
simulation = internal_tide(discretization("AB2-SE"); run_simulation = false)
simulation.Δt  # 196.3 s
run!(simulation)

# Timestep refinement, as in the campaign
internal_tide(discretization("SRK3-SE"); timestep_factor = 1//4, label = "SRK3-SE-dt4")
```

Every simulation carries `VarianceDissipation` callbacks for `b` and `c` and `compute_tracer_dissipation!`, each called
every iteration, and prints progress every 200 iterations. The helpers below are not exported and are reached as
`TimestepperTestCases.<name>`:

```julia
internal_tide_parameters()
internal_tide_grid(; mutable = true)
internal_tide_timestep(Val(timestepper))
internal_tide_stability_parameters()
internal_tide_substeps(barotropic_scheme, scheme = :SplitRungeKutta3; averaging_kernel = OptimizedAsymmetricAveragingKernel())
```

### Campaign

From the repository root:

```sh
julia --project=examples examples/simulations/internal_tide.jl
```

The driver runs the seven production configurations of `discretizations()`, then the timestep refinement of AB2-SE,
WRK3-SE, SRK3-SE and MRK4-SE with `timestep_factor = 1//4` and `1//8` and labels `<label>-dt4` and `<label>-dt8`, for
15 runs in total. A refinement run is skipped when `internal_tide/internal_tide_<TimeStepper>_<label>-dt<n>.jld2`
already exists relative to the working directory. The output schedule is `TimeInterval(1hour)` for every case except
AB2-SE, which saves every `round(Int, 1hour / Δt)` iterations (18, 73 and 147 steps at `Δt`, `Δt / 4` and `Δt / 8`):
an output aligned in time shortens the step that lands on it, and the quasi-AB2 scheme restarts a shortened step with
forward Euler.

`examples/simulations/internal_tide_output_schedule.jl` runs AB2-SE twice, once with `IterationInterval(18)` and once
with `TimeInterval(1hours)`, in `internal_tide/schedule_experiment/{iteration,time}` or under the directory named by
the `AB2_SCHEDULE_OUTPUT` environment variable.
At every step it records the time, the step length, whether the step length changed, and the x-summed profiles of
`Abx`, `Abz` and `Gbz`, and saves them to `every_step.jld2`.

## Output

Each run writes `internal_tide_<timestepper>_<name>.jld2` in the working directory with `JLD2Writer`, overwriting an
existing file. `<timestepper>` is the symbol, for example `SplitRungeKutta3`. `<name>` is the `label` for the
`Discretization` method; for the symbol method it is `split_free_surface` or `implicit_free_surface`, followed by
`_<scheme type>` (for example `_UpwindBiased`) when `tracer_advection` is supplied and the default name is kept.
The published campaign ([Reproducing the paper](../data_examples.md)) holds the 15 files of the driver, such as
`internal_tide_QuasiAdamsBashforth2_AB2-SE.jld2` and `internal_tide_SplitRungeKutta3_WRK3-SE-dt8.jld2`.

| Output | Content | Location |
|:--|:--|:--|
| `u`, `w` | velocity times cell volume | `(Face, Center, Center)`, `(Center, Center, Face)` |
| `b`, `c` | tracer times cell volume | `(Center, Center, Center)` |
| `η` | free-surface displacement | surface |
| `Gbx`, `Gcx` | `(∂x b)²`, `(∂x c)²` times cell volume | `(Face, Center, Center)` |
| `Gbz`, `Gcz` | `(∂z b)²`, `(∂z c)²` times cell volume | `(Center, Center, Face)` |
| `Abx`, `Abz`, `Acx`, `Acz` | advective variance dissipation of `b` and `c` in x and z, from `VarianceDissipation` | x and z faces |
| `Δtb²`, `Δtc²` | `(V b² - (V b²)⁻) / Δt` over the last step, from `compute_tracer_dissipation!` | `(Center, Center, Center)` |

Divide `u`, `w`, `b`, `c` and the squared gradients by the cell volume at their location to recover the fields. Since `c`
starts at zero, its outputs are zero.

## Loading

```julia
load_internal_tide(folder, timestepper, free_surface)
```

reads `joinpath(folder, "internal_tide_$(timestepper)_$(free_surface).jld2")` and returns a `Dict`. It is not exported.

```julia
using TimestepperTestCases
using Oceananigans

d       = discretization("WRK3-SE")
case    = TimestepperTestCases.load_internal_tide("internal_tide", string(d.timestepper), d.label)
refined = TimestepperTestCases.load_internal_tide("internal_tide", string(d.timestepper), d.label * "-dt4")

u = case[:u][end] / case[:VFCC][end]
```

| Key | Content |
|:--|:--|
| `:u`, `:w`, `:b`, `:η` | `FieldTimeSeries` read on disk (stored, volume-weighted values) |
| `:VCCC`, `:VFCC`, `:VCCF` | cell volumes at each output time, recomputed from `η` |
| `:Abx`, `:Abz`, `:Gbx`, `:Gbz`, `:Δtb²` | `FieldTimeSeries` loaded in memory |
| `:abx`, `:abz`, `:abt` | volume-averaged advective dissipation in x, z and their sum, per output time |
| `:gbx`, `:gbz`, `:gbt` | volume-averaged squared buoyancy gradients |
| `:∫Pb`, `:∫Δb²` | domain sums of `Abx + Abz` and of `Δtb²` |
| `:κb` | `FieldTimeSeries` of the local numerical diffusivity `-(Abx + Abz) / 2Gbz` at cell centers |
| `:KE`, `:MKE` | volume-averaged kinetic energy of `(u, w)` and of their x-averages; `KE - MKE` is the TKE in the figures |
| `:η2` | mean of `η²` |
| `:RPE`, `:APE` | reference and available potential energy |

## Results

In the manuscript, AB2-SE is about twice as diffusive as WRK3-SE over the depth, retains less kinetic energy, and
increases the reference potential energy about 1.7 times as much as the Runge-Kutta cases, while WRK3-SE-SM05, SRK3-SE and
MRK4-SE are indistinguishable from WRK3-SE. Under timestep refinement the Runge-Kutta solutions do not change, and
AB2-SE approaches them without matching them at `Δt / 8`. The [manuscript figures](../generated/internal_tide_figures.md)
page reproduces these comparisons from the 15 published files.
