# Diagnostics

This package provides comprehensive diagnostics for quantifying numerical mixing and assessing the performance of time-stepping schemes. The diagnostics are based on exact variance budget methodology that is exact locally in both space and time.

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

The effective numerical diffusivity is computed as:

```math
\kappa_{\text{num}} = \frac{1}{2} \frac{\overline{P_x} + \overline{P_y} + \overline{P_z}}{\overline{|\nabla b|^2}}
```

where the overbar indicates an averaging operator (varies by test case) and `P_x`, `P_y`, `P_z` are the dissipation rates in each direction.

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
case = load_internal_tide("output/", "SplitRungeKutta3", "split_free_surface_DFB")

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

## Interpretation

- **Low RPE increase**: Indicates less irreversible mixing
- **High APE**: Indicates more energy available for instabilities
- **Low numerical diffusivity**: Indicates less spurious mixing from discretization
- **κ_phys / κ_num > 1**: Physical processes dominate over numerical artifacts

As shown in the paper, RK3 schemes consistently show:
- Lower numerical diffusivity than AB2
- Less RPE increase (less irreversible mixing)
- Higher APE (more energy available for instabilities)
- Better preservation of stratification and mesoscale variability

