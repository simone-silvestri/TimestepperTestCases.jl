# Re-entrant channel: configuration and run

The re-entrant channel is the second idealized case of the manuscript. A zonally periodic channel is forced by a wind
stress and a surface buoyancy flux and is restored to a fixed stratification near its northern wall, following the
configurations of Abernathey et al. (2011) and Hill et al. (2012). After a 40-year integration the flow is statistically
equilibrated, so the case measures how each time discretization changes the equilibrated state and how its numerical
mixing compares with the vertical diffusivity of the CATKE parameterization. The figures drawn from the published
outputs are on the [manuscript figures](../generated/channel.md) page.

## Configuration

The values are set in `src/channel_simulation.jl`, with the advection schemes in `src/advection_schemes.jl`.

| Component | Setting |
|:----------|:--------|
| Domain | ``L_x \times L_y \times L_z`` = 1000 km × 2000 km × 3 km, periodic in ``x``, walls at ``y = 0`` and ``y = L_y``, flat bottom |
| Horizontal grid | 200 × 400 cells, ``\Delta x = \Delta y = 5`` km |
| Vertical grid | 90 levels from the surface down: 6 of 10 m, 12 stretched to 39.75 m, 56 of 42 m, 12 stretched back, 4 of 10 m; z-star coordinate |
| Rotation | β-plane with ``f_0 = -10^{-4}`` s⁻¹ and ``\beta = 10^{-11}`` m⁻¹ s⁻¹ |
| Buoyancy | prognostic tracer ``b``, with ``\alpha = 2 \times 10^{-4}`` K⁻¹, ``g = 9.8061`` m s⁻², ``\rho_0 = 999.8`` kg m⁻³, ``c_p = 3994`` J kg⁻¹ K⁻¹ |
| Wind stress | eastward ``\tau_0 \sin(\pi y / L_y)`` with ``\tau_0 = 0.1`` N m⁻², imposed as the top flux ``-(\tau_0 / \rho_0) \sin(\pi y / L_y)`` of ``u`` |
| Surface buoyancy flux | ``Q_0 \cos(3 \pi y / L_y)`` for ``y \le 5 L_y / 6`` and zero further north, positive upward; ``Q_0 = \alpha g \, (10\ \mathrm{W\,m^{-2}}) / (\rho_0 c_p) \approx 4.9 \times 10^{-9}`` m² s⁻³ |
| Sponge restoring | over the northern ``L_y / 10`` (200 km), relaxation time 7 days, toward ``b_\star(z)`` (see below) |
| Bottom drag | quadratic, ``-C_d \, u \lvert \mathbf{u} \rvert`` and ``-C_d \, v \lvert \mathbf{u} \rvert`` in the bottom cell, ``C_d = 1.1 \times 10^{-3}`` |
| Vertical mixing | CATKE (`CATKEEquation(Cᵂϵ = 1.0)`, `minimum_tke = 1e-7`) and a constant viscosity ``\nu = 10^{-4}`` m² s⁻¹; no background diffusivity |
| Horizontal mixing | none beyond the dissipation of the advection schemes |
| Momentum advection | WENO vector invariant: ninth-order vorticity, fifth-order divergence, kinetic-energy gradient and vertical advection |
| Tracer advection | seventh-order WENO, stretched in ``z``; third-order upwind for WRK3-UP |
| Initial condition | equilibrated MITgcm temperature read from `tIni_80y_90L.bin`, converted to ``b = \alpha g T`` |
| Spin-up | 100 days at ``\Delta t = 1`` min with the discretization of the case, after which the clock is reset to zero |
| Integration | 14400 days (40 years of 360 days) at the ``\Delta t`` of the discretization |
| Averaging | consecutive 1800-day (5-year) windows; the last one covers years 35–40 |

The sponge adds to the buoyancy tendency the relaxation term

```math
F_b = \frac{r(y)}{\lambda} \left[ b_\star(z) - b(x, y_N, z) \right] , \qquad
r(y) = \max\left(0, \frac{y - L_y + L_s}{L_s}\right) , \qquad
b_\star(z) = \Delta B \, \frac{e^{z/h} - e^{-L_z/h}}{1 - e^{-L_z/h}} ,
```

with ``L_s = L_y / 10``, ``\lambda = 7`` days, ``\Delta B = 8 \alpha g``, ``h = 1000`` m and ``L_z = 3000`` m. The ramp
``r`` grows linearly from zero at ``y = 0.9 L_y`` to one at the wall. `buoyancy_relaxation` evaluates the anomaly with the
buoyancy of the northernmost row of cells, ``y_N`` (index `j = Ny`), at every latitude of the sponge.

## Discretizations

The seven configurations of the manuscript's discretization table are returned by [`discretizations`](@ref), and a
single one by [`discretization`](@ref) with its label. The label names the output files of a run.

| Label | `timestepper` | Free surface | Substeps, filter | Slow forcing | Tracer advection | ``\Delta t`` |
|:------|:--------------|:-------------|:-----------------|:-------------|:-----------------|:-------------|
| AB2-SE | `:QuasiAdamsBashforth2` | split-explicit | FB, F1 | frozen | WENO, 7th order | 278.6 s (4.6 min) |
| WRK3-SE-SM05 | `:SplitRungeKutta3` | split-explicit | FB, F1 | frozen | WENO, 7th order | 960.2 s (16.0 min) |
| WRK3-SE | `:SplitRungeKutta3` | split-explicit | RK3, F3 | progressive | WENO, 7th order | 960.2 s (16.0 min) |
| SRK3-SE | `:SSPRungeKutta3` | split-explicit | RK3, F3 | frozen | WENO, 7th order | 960.2 s (16.0 min) |
| MRK4-SE | `:ModifiedRungeKutta4` | split-explicit | RK3, F3 | progressive | WENO, 7th order | 1545.0 s (25.75 min) |
| WRK3-IM | `:SplitRungeKutta3` | implicit | none | none | WENO, 7th order | 960.2 s (16.0 min) |
| WRK3-UP | `:SplitRungeKutta3` | split-explicit | RK3, F3 | progressive | upwind, 3rd order | 960.2 s (16.0 min) |

FB and RK3 are the forward–backward and three-stage Runge–Kutta substep integrators. F1 is the filter of Shchepetkin and
McWilliams (`averaging_shape_function`) and F3 is `OptimizedAsymmetricAveragingKernel`. The manuscript and the figures
label WRK3-SE-SM05 as WRK3-SE-F1. WRK3-IM uses `ImplicitFreeSurface()` with its default solver.

The time step follows from [`baroclinic_timestep`](@ref): 0.7 times the stability limit of the first baroclinic mode at
the grid wavenumber ``k = \pi / \Delta x``, that is ``\Delta t = 0.7 \, \theta_\star / (c_1 k)``, with ``\theta_\star``
the imaginary-axis limit of the scheme. The phase speed ``c_1 = \pi^{-1} \int N \, \mathrm{d}z \approx 2.01`` m s⁻¹ is
computed from the stratification of ``b_\star``. Relative to the three-stage schemes (``\theta_\star = \sqrt{3}``), the
ratio of the limits is 0.290 for AB2 and 1.609 for MRK4 with ``a = 0.23``; the environment variable `MRK4_A` sets ``a``
for both the time step and the MRK4 stepper. The barotropic substep count is the smallest multiple of 8 that keeps the
substep Courant number ``c_0 k \Delta\tau`` at or below 0.8 for FB and 1.3 for RK3 substeps, with
``c_0 = \sqrt{g L_z}`` ([`barotropic_substeps`](@ref)). `channel_simulation` logs the time step when it starts.

## Running

[`channel_simulation`](@ref) builds the model, runs the spin-up and the 40-year integration, and returns the
`Simulation` when it has finished. It has a keyword method and a method that takes a `Discretization`:

```julia
channel_simulation(; momentum_advection     = nothing,
                     tracer_advection       = nothing,
                     closure                = default_closure(),
                     zstar                  = true,
                     restart_file           = nothing,
                     free_surface           = nothing,
                     arch                   = GPU(),
                     bottom_height          = nothing,
                     timestepper            = :SplitRungeKutta3,
                     grid                   = default_grid(arch, zstar, bottom_height),
                     initial_file           = "tIni_80y_90L.bin",
                     testcase               = "0",
                     averaging_kernel       = OptimizedAsymmetricAveragingKernel(),
                     barotropic_timestepper = ForwardBackwardScheme(),
                     slow_forcing           = FrozenSlowForcing())

channel_simulation(d::Discretization; kw...)
```

The second method sets `timestepper`, `barotropic_timestepper`, `averaging_kernel`, `slow_forcing`, `free_surface`,
`tracer_advection` and `testcase = d.label` from `d`, and passes the remaining keywords through. It is the method that
reproduces a row of the table above; the keyword defaults pair forward–backward substeps with the F3 filter and frozen
forcing. With `free_surface = nothing` a split-explicit free surface is built with the substep count described above.
With `bottom_height`, a function `(x, y) -> z`, the grid becomes an `ImmersedBoundaryGrid` and the quadratic drag also
acts on the immersed bottom. The campaign uses the flat bottom.

The default `initial_file`, a path relative to the working directory, holds a big-endian `Float64` array of
200 × 400 × 90 temperatures. The file is distributed neither with the package nor with the Zenodo record;
`initial_file = nothing` starts from ``b_\star(z)`` instead. The default architecture is `GPU()`, which needs CUDA loaded:

```julia
using CUDA
using TimestepperTestCases

simulation = channel_simulation(discretization("WRK3-SE"); initial_file = nothing)
```

Without a `restart_file`, the run starts with the 100-day spin-up, which ends with `restart<label>_iteration144000.jld2`.
With `restart_file`, the state is loaded from that checkpoint, the spin-up is skipped, and the clock is reset to zero
before the 40-year integration. A restart file therefore sets the initial state of the integration, and a run started
from a later `channel_checkpoint_<label>` file integrates another 14400 days from that state.

```julia
simulation = channel_simulation(discretization("MRK4-SE");
                                restart_file = "restartMRK4-SE_iteration144000.jld2")
```

`examples/simulations/channel.jl` runs the campaign on the GPU from the repository root, one configuration after the
other, writing to the working directory. Without arguments it runs all seven labels in the order of the table;
otherwise it runs the labels given. If `CHANNEL_RESTART_DIR` is set, each case starts from
`restart<label>_iteration144000.jld2` in that directory.

```sh
julia --project=. examples/simulations/channel.jl
julia --project=. examples/simulations/channel.jl WRK3-SE MRK4-SE
CHANNEL_RESTART_DIR=/path/to/restarts julia --project=. examples/simulations/channel.jl
```

## Output

Each run writes, with `<label>` the `testcase` string:

| File | Content | Schedule |
|:-----|:--------|:---------|
| `snapshots_<label>.jld2` | `u`, `v`, `w`, `b`, `e`, `Abx`, `Aby`, `Abz`, `Dbx`, `Dby`, `Dbz`, `Gbx`, `Gby`, `Gbz`, `η`, `VCCC`, `VFCC`, `VCFC`, `VCCF` | every 360 days, at two consecutive iterations |
| `averages_<label>.jld2` | the snapshot fields and `vb`, `wb`, `u²`, `v²` | averages over consecutive 1800-day windows |
| `mixing_<label>.jld2` | `b`, `η` | every 360 days |
| `channel_checkpoint_<label>_iteration<N>.jld2` | model state | every 1800 days |
| `restart<label>_iteration<N>.jld2` | model state, `N = 144000` at the end of the spin-up | every 100 days of the spin-up |

`e` is the CATKE turbulent kinetic energy. `Abx`, `Aby` and `Abz` are the advective and `Dbx`, `Dby` and `Dbz` the
diffusive contributions to the buoyancy-variance budget in each direction, computed at every iteration by Oceananigans'
`VarianceDissipation`. `Gbx`, `Gby` and `Gbz` are the squared buoyancy gradients multiplied by the volume of the cell on
which they are located, and `VCCC`, `VFCC`, `VCFC` and `VCCF` are the cell volumes at the four staggered locations,
which vary in time under z-star. The ratios of these terms give the numerical and physical diffusivities described in
[diagnostics](../diagnostics.md). The mixing file keeps the instantaneous `b` and `η` because the reference potential
energy is not linear in ``b`` and cannot be computed from the averages.

Zenodo record 24315107 ([doi:10.5281/zenodo.24315107](https://doi.org/10.5281/zenodo.24315107)) holds
`averages_<label>.jld2` and `mixing_<label>.jld2` for the seven labels. Snapshots, checkpoints and restart files are
not published. Downloading the files and running the figure example are described in
[Reproducing the paper](../data_examples.md).

## Loading

`TimestepperTestCases.load_channel(folder, case_number; arch = CPU())` reads `snapshots_<case_number>.jld2` from
`folder`, so it applies to local runs. It returns a `Dict` with `:u`, `:v`, `:w`, `:b` and `:η` as on-disk
`FieldTimeSeries`, the cell volumes `:VCCC`, `:VFCC`, `:VCFC` and `:VCCF` recomputed from `η`, and one value per
snapshot of the volume-mean kinetic energy `:KE`, the kinetic energy of the zonal mean `:MKE`, the mean of ``\eta^2``
`:η2`, and the reference and available potential energies `:RPE` and `:APE`.

```julia
using Oceananigans
using TimestepperTestCases

case = TimestepperTestCases.load_channel("path/to/run", "WRK3-SE"; arch = CPU())
case[:KE]
```

The averages and mixing files are Oceananigans JLD2 output and are read with `FieldTimeSeries`. The figure example
takes zonal means of the last averaging window through `prepare_channel` in `examples/shared/compact_data.jl`.

```julia
using Oceananigans

Abx = FieldTimeSeries("averages_WRK3-SE.jld2", "Abx"; backend = OnDisk())
```

## Results

Averaged over years 35–40, the horizontal numerical diffusivity is of order 1–10 m² s⁻¹ in every configuration and
the total numerical diffusivity ranges from ``10^{-6}`` to ``10^{-4}`` m² s⁻¹. WRK3-UP mixes the most, with
``\kappa_{\rm phys} / \kappa_{\rm num} \lesssim 1`` over the whole water column. AB2-SE mixes less than WRK3-UP and more
than WRK3-SE and the other Runge–Kutta configurations with WENO tracer advection. In every configuration other than
WRK3-UP the ratio exceeds one, with local peaks above 10 in the interior, and it is largest for the Runge–Kutta
schemes. The horizontal contribution dominates the variance destruction, while the vertical contribution has patches
of variance production. The profiles and sections are on the [manuscript figures](../generated/channel.md) page.
