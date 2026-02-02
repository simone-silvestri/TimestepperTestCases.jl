# Linear Stability Analysis

This package includes tools for analyzing the linear stability of coupled barotropic-baroclinic systems with different time-stepping schemes. The stability analysis extends the modal decomposition approach to multi-stage Runge-Kutta schemes.

## Overview

The stability analysis decomposes the system into barotropic and baroclinic modes and computes eigenvalues of the evolution matrix. For a two-dimensional (x, z) hydrostatic system with no rotation, the analysis considers:

- Barotropic mode: Depth-independent mode
- First baroclinic mode: Projection onto vertical mode `φ = √2 cos(πz/H)`

## Key Findings

The analysis reveals critical constraints for split-explicit RK implementations:

1. **Transport Velocity Choice**: Split-explicit RK schemes require the averaged barotropic transport velocity (`w★`) for tracer advection, not the previous substep's velocity (`w^m`), to maintain baroclinic stability.

2. **Stability Comparison**: RK3-SE with correct coupling shows baroclinic mode amplification similar to the single-mode RK analysis, while incorrect coupling (`w^m`) leads to instabilities.

3. **Implicit vs. Split-Explicit**: RK3-IM shows increased dissipation compared to RK3-SE for both barotropic and baroclinic modes, as expected from the implicit free surface treatment.

## Using the Stability Analysis

The stability analysis is implemented in the `linear_stability.ipynb` notebook. Key functions compute:

- Evolution matrices for different timesteppers
- Eigenvalues for barotropic and baroclinic modes
- Amplification factors as functions of CFL number and stratification

### Running the Analysis

```julia
# See linear_stability.ipynb for full implementation
# The notebook computes eigenvalues for:
# - AB2-SE
# - RK3-SE with w^m (incorrect)
# - RK3-SE with w★ (correct)
# - RK3-IM
```

## Mathematical Framework

The analysis linearizes the primitive equations around a rest state and decomposes variables as:

```math
Y(x, z, t) = \sum_{q=0}^\infty Y_q(x, t) M_q(z)
```

where `M_0 = 1` (barotropic) and `M_1 = √2 cos(πz/H)` (first baroclinic).

The discrete evolution matrix is constructed from the time-stepping scheme, and eigenvalues determine stability:

- `|λ| < 1`: Stable (dissipative)
- `|λ| = 1`: Neutral (non-dissipative)
- `|λ| > 1`: Unstable

## Results

The stability analysis shows:

1. **Baroclinic Mode**: RK3-SE with `w★` maintains stability similar to single-mode RK, while `w^m` leads to growing instabilities with increasing stratification.

2. **Barotropic Mode**: Both RK3-SE and RK3-IM maintain stable barotropic modes, with RK3-IM showing more dissipation.

3. **CFL Dependence**: RK3 schemes remain stable over larger CFL ranges than AB2, consistent with the amplification factor analysis.

## References

The analysis methodology is described in detail in the paper, extending the approach of Demange et al. (2019) to multi-stage schemes. See the paper's Section 4 for the complete derivation.

