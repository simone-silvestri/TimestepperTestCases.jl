# Getting Started

## Installation

Install TimestepperTestCases.jl using Julia's package manager:

```julia
using Pkg
Pkg.add("TimestepperTestCases")
```

Or add it directly from the GitHub repository:

```julia
using Pkg
Pkg.add(url="https://github.com/simone-silvestri/TimestepperTestCases.jl")
```

## Basic Usage

After installation, load the package:

```julia
using TimestepperTestCases
```

### Running a Test Case

The package provides three main test case functions:

#### Internal Tide

```julia
using TimestepperTestCases

# Run with RK3 and split-explicit free surface
sim = internal_tide(:SplitRungeKutta3)

# Or with AB2
sim = internal_tide(:QuasiAdamsBashforth2)
```

#### Channel Simulation

```julia
# Run channel simulation with default settings
sim = channel_simulation(timestepper=:SplitRungeKutta3)

# Or customize closure and grid
sim = channel_simulation(
    timestepper=:SplitRungeKutta3,
    closure=default_closure(),
    zstar=true
)
```

### Loading Simulation Output

After running a simulation, load the output data:

```julia
# For internal tide
case = load_internal_tide("output_folder/", "SplitRungeKutta3", "split_free_surface_DFB")

# Access diagnostics
case[:KE]   # Kinetic energy time series
case[:RPE]  # Reference potential energy
case[:APE]  # Available potential energy
case[:κb]   # Numerical diffusivity
```

## Dependencies

TimestepperTestCases.jl requires:

- Julia 1.9 or later
- Oceananigans.jl (for the ocean model)
- JLD2 (for reading/writing simulation output)
- Various plotting and analysis packages (see `Project.toml`)

## Next Steps

- Learn about the [experiments](experiments/internal_tide.md) in detail
- Understand the [diagnostics](diagnostics.md) used to quantify numerical mixing
- Run the [Julia examples](data_examples.md) with local or published data

