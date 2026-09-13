#####
##### The global tripolar configuration: the near-global machinery of `near_global.jl` on the NEMO eORCA025
##### mesh, which carries its own metrics and its own bathymetry.
#####

using NumericalEarth

"""
    global_ocean_grid(arch = CPU(); Nz, depth, surface_spacing, z, dataset, major_basins, halo)

The eORCA025 tripolar grid, on the vertical coordinate of `near_global_vertical_discretization`. The mesh
carries its own metrics and its own bathymetry; the northern boundary is the `RightFaceFolded` seam and the
zonal direction is periodic.

$(SIGNATURES)
"""
function global_ocean_grid(arch = CPU();
                           Nz = 100,
                           depth = 5000meters,
                           surface_spacing = 5meters,
                           z = near_global_vertical_discretization(Nz, depth, surface_spacing),
                           dataset = ORCAQuarter(),
                           major_basins = 1,
                           halo = (7, 7, 7))

    return ORCAGrid(arch, Oceananigans.defaults.FloatType;
                    dataset, z, Nz, halo, major_basins, active_cells_map = true)
end

"""
    global_ocean(discretization = :SplitRungeKutta3; arch, grid, kw...)

Run `near_global` on the tripolar grid of `global_ocean_grid`, under the `global_ocean` output prefix.
`discretization` is a timestepper symbol or a `Discretization`, and every keyword of `near_global` applies.
"""
global_ocean(discretization = :SplitRungeKutta3; arch = CPU(),
             grid = global_ocean_grid(arch), kw...) = near_global(discretization; arch, grid, prefix = "global_ocean", kw...)
