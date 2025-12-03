# Internal Tide Test Case

The internal tide test case simulates tidal flow over a Gaussian seamount in a stratified ocean. This configuration isolates the role of time discretization, as spatial advection plays a secondary role in this mostly linear regime.

## Configuration

The domain spans 2000 km horizontally (periodic) and 2 km vertically (bounded), with a Gaussian seamount of height 250 m and width 20 km. The domain is initially stratified with `N² = 10⁻⁴ s⁻²` and forced by an oscillatory tidal forcing.

### Parameters

```julia
using TimestepperTestCases

# Get default parameters
params = internal_tide_parameters()
```

Key parameters:
- Domain: 2000 km × 2 km
- Grid: 256 × 128 points
- Seamount: 250 m height, 20 km width
- Tidal period: 12.421 hours (M₂ tide)
- Initial stratification: `N² = 10⁻⁴ s⁻²`
- Time step: 5 minutes (AB2) or 10 minutes (RK3)
- Simulation duration: 40 days

## Running the Simulation

```julia
using TimestepperTestCases

# Run with RK3 and split-explicit free surface
sim = internal_tide(:SplitRungeKutta3)

# Or with implicit free surface
sim = internal_tide(:SplitRungeKutta3; 
                    free_surface=ImplicitFreeSurface(internal_tide_grid()))

# Run with AB2 for comparison
sim = internal_tide(:QuasiAdamsBashforth2)
```

## Output Fields

The simulation outputs:
- `u`, `w`: Horizontal and vertical velocities
- `b`, `c`: Buoyancy and passive tracer
- `η`: Free surface elevation
- `Δtc²`, `Δtb²`: Tracer variance dissipation rates
- `Gbx`, `Gbz`, `Gcx`, `Gcz`: Gradient squared fields

## Expected Results

As shown in the paper, RK3 schemes retain more energy than AB2:
- Sharper velocity structures (`u′`, `w`)
- Better preserved stratification (`N²`)
- Lower numerical diffusivity throughout the water column
- Higher kinetic energy and available potential energy

## Loading and Analyzing Results

```julia
# Load results
case = load_internal_tide("output_folder/", "SplitRungeKutta3", "split_free_surface_DFB")

# Access time series
case[:KE]   # Kinetic energy
case[:RPE]  # Reference potential energy (irreversible mixing)
case[:APE]  # Available potential energy
case[:abx]  # Volume-averaged x-direction dissipation
case[:abz]  # Volume-averaged z-direction dissipation
case[:κb]   # Effective numerical diffusivity
```

## Physical Interpretation

This test case demonstrates that temporal discretization alone can introduce significant numerical mixing, with AB2 showing diffusivity comparable to low-order spatial schemes. The RK3 formulation substantially reduces this mixing, preserving wave energy and stratification over the 40-day integration.

