# Idealized Coastal Baroclinic Adjustment

The idealized coast test case simulates baroclinic adjustment in a coastal domain with a sloping bottom and initial meridional salinity gradient. This configuration introduces nonlinearities and demonstrates how numerical mixing can suppress submesoscale variability.

## Configuration

The domain is 192 km × 192 km × 103 m with a linearly sloping bathymetry (5 m deep at the southern edge, reaching full depth at 96 km). The initial hydrography features a meridional salinity gradient that drives baroclinic instabilities.

### Parameters

```julia
using TimestepperTestCases

# Run with default settings (forced, high resolution)
sim = idealized_coast(:SplitRungeKutta3)

# Low resolution version
sim = idealized_coast(:SplitRungeKutta3; lowres=true)

# Unforced version (no wind stress)
sim = idealized_coast(:SplitRungeKutta3; forced=false)
```

Key parameters:
- Domain: 192 km × 192 km × 103 m
- Grid: 250 × 250 × 40 points (high res) or 96 × 96 × 40 (low res)
- Bottom slope: Linear from 5 m to 103 m depth
- Initial stratification: `N² = 10⁻⁴ s⁻²`
- Meridional gradient: `M² = 1.2×10⁻⁶ s⁻²` (for y < 50 km)
- Time step: 5 minutes (AB2) or 10 minutes (RK3)
- Simulation duration: 40 days
- Wind forcing: Oscillatory stress after 4-day spin-up (if `forced=true`)

## Running the Simulation

```julia
using TimestepperTestCases

# Standard forced case
sim = idealized_coast(:SplitRungeKutta3; forced=true, lowres=false)

# Customize free surface
sim = idealized_coast(:SplitRungeKutta3; 
                      free_surface=ImplicitFreeSurface(grid))
```

## Output Fields

The simulation outputs:
- `u`, `v`, `w`: Velocity components
- `T`, `S`: Temperature and salinity
- `b`: Buoyancy
- `η`: Free surface elevation
- `Abx`, `Aby`, `Abz`: Advective buoyancy dissipation (x, y, z)
- `ATx`, `ATy`, `ATz`: Advective temperature dissipation
- `ASx`, `ASy`, `ASz`: Advective salinity dissipation
- `Gbx`, `Gby`, `Gbz`: Buoyancy gradient squared
- `GTx`, `GTy`, `GTz`: Temperature gradient squared
- `GSx`, `GSy`, `GSz`: Salinity gradient squared
- `κu`, `κc`: CATKE diffusivities (if forced)

## Expected Results

As shown in the paper:
- RK3-SE maintains sharper stratification and stronger submesoscale variability
- AB2-SE shows reduced eddy activity and weaker stratification
- Horizontal dissipation dominates over vertical dissipation
- RK3-SE preserves more kinetic energy variance across scales

## Loading and Analyzing Results

```julia
# Load results
case = load_idealized_coast("output_folder/", "CATKE", "split_free_surface_", "SplitRungeKutta3")

# Access diagnostics
case[:KE]   # Kinetic energy
case[:MKE]  # Mean kinetic energy
case[:RPE]  # Reference potential energy
case[:APE]  # Available potential energy
case[:abx], case[:aby], case[:abz]  # Volume-averaged dissipation
case[:gbx], case[:gby], case[:gbz]  # Volume-averaged gradient squared
```

## Physical Interpretation

This test case demonstrates that numerical mixing can suppress physical instabilities. The reduced dissipation in RK3-SE allows baroclinic instabilities to develop more fully, preserving submesoscale variability that AB2-SE suppresses. The horizontal dissipation dominates because eddy activity occurs primarily at the surface where horizontal gradients are strongest.

