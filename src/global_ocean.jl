#####
##### The global tripolar configuration: the near-global machinery of `near_global.jl` on the NEMO eORCA025
##### mesh, which carries its own metrics and its own bathymetry.
#####

using NumericalEarth

"""
    global_ocean_grid(arch = CPU(); Nz, depth, surface_spacing, z, dataset, major_basins, minimum_depth, halo)

The eORCA025 tripolar grid, on the vertical coordinate of `near_global_vertical_discretization`. The mesh
carries its own metrics and its own bathymetry; the northern boundary is the `RightFaceFolded` seam and the
zonal direction is periodic. Columns shallower than `minimum_depth` are land, as in `near_global_grid`.

$(SIGNATURES)
"""
function global_ocean_grid(arch = CPU();
                           Nz = 100,
                           depth = 5000meters,
                           surface_spacing = 5meters,
                           z = near_global_vertical_discretization(Nz, depth, surface_spacing),
                           dataset = ORCAQuarter(),
                           major_basins = 1,
                           minimum_depth = 30meters,
                           halo = (7, 7, 7))

    grid = ORCAGrid(arch, Oceananigans.defaults.FloatType;
                    dataset, z, Nz, halo, major_basins, active_cells_map = true)

    fill_degenerate_metrics!(grid)

    bottom_height = grid.immersed_boundary.bottom_height
    parent(bottom_height) .= land_above_minimum_depth.(parent(bottom_height), minimum_depth)

    return ImmersedBoundaryGrid(grid.underlying_grid, GridFittedBottom(bottom_height); active_cells_map = true)
end

@inline land_above_minimum_depth(h, minimum_depth) = ifelse(h > -minimum_depth, zero(h), h)

const horizontal_metric_names = (:Δxᶜᶜᵃ, :Δxᶠᶜᵃ, :Δxᶜᶠᵃ, :Δxᶠᶠᵃ,
                                 :Δyᶜᶜᵃ, :Δyᶠᶜᵃ, :Δyᶜᶠᵃ, :Δyᶠᶠᵃ,
                                 :Azᶜᶜᵃ, :Azᶠᶜᵃ, :Azᶜᶠᵃ, :Azᶠᶠᵃ)

"""
    fill_degenerate_metrics!(grid)

Replace every nonpositive horizontal spacing and area of `grid` with the nearest positive value on the same row
(or, for a row with none, the nearest row that has one). The eORCA025 mesh carries zero metrics on the rows
under the northern fold and in the southern halo, whose inverses put `0 × ∞` into the continuity and momentum
operators.
"""
function fill_degenerate_metrics!(grid)
    underlying_grid = grid isa Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid ? grid.underlying_grid : grid

    for name in horizontal_metric_names
        metric = parent(getproperty(underlying_grid, name))
        host = Array(metric)
        filled = fill_degenerate_rows!(host)
        filled > 0 && @info "fill_degenerate_metrics!: $filled nonpositive entries of $name replaced"
        copyto!(metric, host)
    end

    return grid
end

function fill_degenerate_rows!(metric)
    filled = 0
    valid_rows = Int[]

    for j in axes(metric, 2)
        row = view(metric, :, j)
        valid = findall(>(0), row)
        isempty(valid) && continue
        push!(valid_rows, j)

        for i in eachindex(row)
            if !(row[i] > 0)
                row[i] = row[valid[argmin(abs.(valid .- i))]]
                filled += 1
            end
        end
    end

    for j in axes(metric, 2)
        j in valid_rows && continue
        nearest = valid_rows[argmin(abs.(valid_rows .- j))]
        filled += count(x -> !(x > 0), view(metric, :, j))
        metric[:, j] .= metric[:, nearest]
    end

    return filled
end

"""
    global_ocean(discretization = :SplitRungeKutta3; arch, grid, kw...)

Run `near_global` on the tripolar grid of `global_ocean_grid`, under the `global_ocean` output prefix.
`discretization` is a timestepper symbol or a `Discretization`, and every keyword of `near_global` applies.
"""
global_ocean(discretization = :SplitRungeKutta3; arch = CPU(), grid = global_ocean_grid(arch),
             prefix = "global_ocean", kw...) = near_global(discretization; arch, grid, prefix, kw...)

"""
    global_ocean_cost(discretization = :SplitRungeKutta3; arch, grid, warmup_steps, measured_steps, kw...)

A short `global_ocean` run that only measures the cost per step of `discretization`. The barotropic substep
count is sized for the nominal time step of the scheme, but the model steps at a third of it, the cold-start step
of `near_global`: the work per step is the same, since neither the substep count, nor CATKE, nor the adaptive
vertically implicit advection depend on the step, and the run stays far from its stability limit.
`warmup_steps` absorb the compilation and are not timed; `seconds_per_step` is measured over `measured_steps`.

The implicit free surface is the exception: the iterations of its solver grow with the step, so its
`measured_steps` run at the nominal time step.

No output is written and the buoyancy-variance budget is off. The record goes to
`global_ocean_cost_<label>_cost.jld2`.
"""
function global_ocean_cost(discretization = :SplitRungeKutta3; arch = CPU(), grid = global_ocean_grid(arch),
                           warmup_steps = 200, measured_steps = 1000, kw...)

    timestepper = discretization isa Discretization ? discretization.timestepper : discretization
    implicit_free_surface = discretization isa Discretization && discretization.implicit_free_surface

    nominal_timestep  = near_global_timestep(Val(timestepper))
    warmup_timestep   = nominal_timestep / 3
    measured_timestep = implicit_free_surface ? nominal_timestep : warmup_timestep

    cold_start_duration = warmup_steps * warmup_timestep
    stop_time = cold_start_duration + measured_steps * measured_timestep

    return global_ocean(discretization; arch, grid, prefix = "global_ocean_cost",
                        Δt = measured_timestep, nominal_timestep, cold_start_Δt = warmup_timestep,
                        cold_start_duration, stop_time, cost_stop_time = stop_time,
                        dissipation = false, output = false, progress_interval = IterationInterval(100), kw...)
end
