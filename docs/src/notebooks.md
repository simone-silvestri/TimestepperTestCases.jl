# Notebooks

This repository includes Jupyter notebooks for reproducing figures and analysis from the paper. All notebooks are located in the `notebooks/` directory.

## Available Notebooks

### Figure Generation

- **`figure_1.ipynb`**: Generates Figure 1 from the paper (amplification factors)
- **`figure_2.ipynb`**: Generates Figure 2 from the paper (linear stability analysis)

### Results Analysis

- **`Results_internal_tide.ipynb`**: Analysis and visualization of internal tide results
- **`Results_idealized_coast.ipynb`**: Analysis and visualization of idealized coast results
- **`Results_channel.ipynb`**: Analysis and visualization of channel simulation results

### Stability Analysis

- **`linear_stability.ipynb`**: Linear stability analysis for coupled barotropic-baroclinic modes

## Running Notebooks

### Prerequisites

Install required packages:

```julia
using Pkg
Pkg.add(["IJulia", "Plots", "CairoMakie", "JLD2"])
```

### Starting Jupyter

```julia
using IJulia
notebook()
```

Or from the command line:

```bash
julia -e "using IJulia; notebook()"
```

### Running a Notebook

1. Navigate to the `notebooks/` directory
2. Open the desired notebook in Jupyter
3. Ensure simulation output files are in the expected locations (see notebook comments)
4. Run all cells

## Notebook Descriptions

### Figure 1: Amplification Factors

This notebook reproduces the amplification factor analysis comparing AB2 and RK schemes for linear advection. It shows how RK schemes maintain amplification factors closer to unity over a wider CFL range, indicating less implicit dissipation.

### Figure 2: Linear Stability

This notebook reproduces the linear stability analysis showing eigenvalues for the coupled barotropic-baroclinic system. It demonstrates the importance of using the averaged barotropic transport velocity (`w★`) for tracer advection in split-explicit RK schemes.

### Results Notebooks

These notebooks load simulation output and generate diagnostic plots including:
- Time series of kinetic energy, RPE, APE
- Vertical profiles of numerical diffusivity
- Spatial snapshots of velocity, buoyancy, and dissipation fields
- Power spectra
- Mean fields and differences between timesteppers

### Linear Stability Notebook

This notebook implements the modal decomposition analysis described in the paper, computing eigenvalues for different coupling strategies and demonstrating stability constraints for multi-stage schemes.

## Data Requirements

Most notebooks require simulation output files. Run the simulations first:

```julia
using TimestepperTestCases

# Run simulations
internal_tide(:SplitRungeKutta3)
internal_tide(:QuasiAdamsBashforth2)

idealized_coast(:SplitRungeKutta3)
idealized_coast(:QuasiAdamsBashforth2)

channel_simulation(timestepper=:SplitRungeKutta3)
channel_simulation(timestepper=:QuasiAdamsBashforth2)
```

Output files should be in the directories expected by the notebooks (typically `internal_tide/`, `idealized_coast/`, and `channel-simulation/`).

## Customization

Notebooks can be customized to:
- Compare different timesteppers
- Analyze different time periods
- Generate custom visualizations
- Compute additional diagnostics

Modify the notebook code as needed for your analysis.

