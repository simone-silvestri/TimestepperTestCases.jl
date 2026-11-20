# Variance-budget diagnostics

Numerical mixing is quantified from the discrete tracer-variance budget, which closes exactly in space and time for each
cell and each timestep. The diagnosed variance destruction, divided by the squared tracer gradient, defines the numerical
diffusivity shown in the manuscript figures.

## Variance Budget Methodology

The variance budget approach quantifies numerical mixing by tracking the dissipation of tracer variance. For a tracer `c`, the variance budget reads:

```math
\frac{\partial c^2}{\partial t} = -2c \nabla \cdot (\mathbf{u} c) + \text{sources} - \text{sinks}
```

The numerical dissipation `P` is computed as:

```math
P = 2 F_\star \delta c - U_\star \delta c^2
```

where `F_\star` is the effective advective flux and `U_\star` is the effective advective velocity, computed differently for AB2 and RK schemes.

### For AB2

```math
F_\star = \left(\frac{3}{2} + \epsilon\right) F^n - \left(\frac{1}{2} + \epsilon\right) F^{n-1}
```

```math
U_\star = \left(\frac{3}{2} + \epsilon\right) U^n - \left(\frac{1}{2} + \epsilon\right) U^{n-1}
```

### For RK Schemes

For low-storage RK schemes, only the last substep contributes to the budget:

```math
F_\star = F^{M-1}, \quad U_\star = U^{M-1}
```

where `M-1` indicates the second-to-last substep.

## Numerical Diffusivity

The numerical diffusivity of a region is

```math
\kappa_{\text{num}} = -\frac{\sum (P_x + P_y + P_z)}{2 \sum |\nabla b|^2},
```

where `P_x`, `P_y`, `P_z` are the variance tendencies in each direction (saved as `Abx`, `Aby`, `Abz`), `|∇b|²` is the
squared buoyancy gradient (saved as `Gbx`, `Gby`, `Gbz`), and both sums run over the cells and output windows of the region.
Both saved quantities already include cell volumes, so the sums are volume integrals. The ratio is formed after summing,
never cell by cell. Variance destruction is negative, so κ > 0 is diffusive and κ < 0 anti-diffusive.

## Reference Potential Energy (RPE)

Reference potential energy represents the minimum potential energy achievable by adiabatic re-sorting of the density field. It quantifies irreversible mixing:

```math
\text{RPE} = \int_V \rho z^\star dV
```

where `z^\star` is the re-sorted height computed by sorting buoyancy and computing cumulative volume distribution.

## Available Potential Energy (APE)

Available potential energy is the difference between actual and reference potential energy:

```math
\text{APE} = \int_V \rho (z - z^\star) dV
```

APE represents the energy available for conversion to kinetic energy through baroclinic instabilities.

## Using Diagnostics

### During Simulation

Add diagnostics as callbacks:

```julia
using TimestepperTestCases
using Oceananigans.Models.VarianceDissipationComputations

# For buoyancy variance dissipation
ϵb = BuoyancyVarianceDissipation(grid)
add_callback!(simulation, ϵb, IterationInterval(1))

# For tracer variance dissipation (internal tide)
add_callback!(simulation, compute_tracer_dissipation!, IterationInterval(1))
```

### Post-Processing

Load and compute diagnostics from saved output:

```julia
# Load case
case = load_internal_tide("internal_tide/", "SplitRungeKutta3", "WRK3-SE")

# Access pre-computed diagnostics
case[:RPE]  # Reference potential energy time series
case[:APE]  # Available potential energy time series
case[:κb]   # Numerical diffusivity FieldTimeSeries
case[:abx]  # Volume-averaged x-direction dissipation
case[:abz]  # Volume-averaged z-direction dissipation
```

### Computing RPE/APE Manually

```julia
using TimestepperTestCases

# From a case dictionary
rpe_ape = compute_rpe_density(case)
rpe_ape.rpe  # RPE time series
rpe_ape.ape  # APE time series

# From fields directly
zstar = calculate_z★_diagnostics(b_field, vol_field)
rpe_ape_fields = compute_rpe_density(b_field, vol_field)
```

## Results

The diagnosed diffusivities of each experiment are shown on the figure pages of the
[internal tide](generated/internal_tide_figures.md), the [re-entrant channel](generated/channel.md) and the
[quarter-degree global ocean](generated/global_ocean_results.md).
