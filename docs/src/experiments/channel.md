# Equilibrated Channel Flow

The channel simulation is an idealized re-entrant periodic channel forced by wind stress and surface heat flux. Unlike the transient configurations, this case relaxes toward statistical equilibrium, allowing assessment of how time-stepping schemes affect the equilibrated solution and how numerical mixing interacts with parameterized physical diffusivity.

## Configuration

The domain is 1000 km × 2000 km × 3 km with periodic zonal boundaries and bounded meridional boundaries. The channel is forced by sinusoidal wind stress and variable surface heat flux, with buoyancy restored at the northern boundary.

### Parameters

```julia
using TimestepperTestCases

# Run with default settings
sim = run_channel_simulation(timestepper=:SplitRungeKutta3)

# Customize closure
sim = run_channel_simulation(
    timestepper=:SplitRungeKutta3,
    closure=default_closure(),
    zstar=true
)

# Use custom grid
sim = run_channel_simulation(
    timestepper=:SplitRungeKutta3,
    grid=custom_grid
)
```

Key parameters:
- Domain: 1000 km × 2000 km × 3 km
- Grid: 200 × 400 × 90 points (non-uniform vertical spacing)
- Time step: 5 minutes (AB2), 10 minutes (RK3), or 20 minutes (RK6)
- Simulation duration: 40 years (last 5 years averaged)
- Closure: CATKE with background diffusivity `κ = 10⁻⁵ m²/s`, `ν = 10⁻⁴ m²/s`
- Wind stress: Sinusoidal profile `τ = -0.1/ρ * sin(π * y / Ly)`
- Heat flux: Variable profile `Qᵇ * cos(3π * y / Ly)` for `y < 5/6 * Ly`

## Running the Simulation

```julia
using TimestepperTestCases

# Standard run
sim = run_channel_simulation(timestepper=:SplitRungeKutta3)

# Restart from checkpoint
sim = run_channel_simulation(
    timestepper=:SplitRungeKutta3,
    restart_file="channel_checkpoint_0.jld2"
)
```

## Output Fields

The simulation outputs:
- `u`, `v`, `w`: Velocity components
- `b`: Buoyancy
- `η`: Free surface elevation
- `Abx`, `Aby`, `Abz`: Advective buoyancy dissipation
- `Gbx`, `Gby`, `Gbz`: Buoyancy gradient squared
- Snapshots: Every 360 days
- Averages: 5-year averages

## Expected Results

As shown in the paper:
- RK3-SE achieves `κ_phys / κ_num > 2` throughout most of the water column
- AB2-SE hovers around `κ_phys / κ_num ≈ 1`
- RK3-UP (with low-order advection) shows `κ_phys / κ_num < 1`
- RK3-SE produces cooler equilibrated solutions (less spurious mixing)
- Horizontal numerical diffusivity dominates over vertical

## Loading and Analyzing Results

```julia
# Load channel results
case = load_channel("output_folder/", "0"; arch=CPU())

# Access diagnostics
case[:KE]   # Kinetic energy time series
case[:MKE]  # Mean kinetic energy
case[:RPE]  # Reference potential energy
case[:APE]  # Available potential energy
case[:η2]   # Free surface variance
```

## Physical Interpretation

This test case demonstrates the critical importance of maintaining `κ_phys / κ_num > 1` for long-term climate integrations. When numerical mixing overwhelms physical mixing (as in RK3-UP or near-surface AB2-SE), the equilibrated solution is significantly warmer and less stratified. RK3-SE maintains physical processes as dominant throughout most of the water column, which is essential for preserving stratification and mesoscale variability in climate applications.

