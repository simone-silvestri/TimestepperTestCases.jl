"""
    load_internal_tide(folder, timestepper, free_surface)

Load and process internal tide simulation output data.

$(SIGNATURES)

# Arguments
- `folder`: Directory containing the simulation output files
- `timestepper`: Timestepper name string (e.g., `"SplitRungeKutta3"`)
- `free_surface`: Free surface type string (e.g., `"split_free_surface_DFB"`)

# Returns
- Dictionary containing:
  - `:u`, `:w`, `:b`, `:η`: Velocity, buoyancy, and free surface FieldTimeSeries
  - `:VCCC`, `:VFCC`, `:VCCF`: Volume FieldTimeSeries at different grid locations
  - `:Abx`, `:Abz`: Advective buoyancy dissipation in x and z directions
  - `:Gbx`, `:Gbz`: Buoyancy gradient squared in x and z directions
  - `:abx`, `:abz`, `:abt`: Volume-averaged dissipation rates
  - `:gbx`, `:gbz`, `:gbt`: Volume-averaged gradient squared
  - `:κb`: Effective numerical diffusivity
  - `:KE`, `:MKE`: Kinetic energy and mean kinetic energy time series
  - `:η2`: Free surface variance time series
  - `:RPE`, `:APE`: Reference and available potential energy time series

This function loads simulation output, computes volume fields accounting for free surface
variations, and calculates various diagnostics used for analyzing numerical mixing.
"""
function load_internal_tide(folder, timestepper, free_surface)
    path = folder * "internal_tide_" * timestepper * "_" * free_surface * ".jld2"
    case = Dict()

    grid = internal_tide_grid()
    case[:u] = FieldTimeSeries(path, "u"; backend=OnDisk())
    case[:w] = FieldTimeSeries(path, "w"; backend=OnDisk())
    case[:b] = FieldTimeSeries(path, "b"; backend=OnDisk())
    case[:η] = FieldTimeSeries(path, "η"; backend=OnDisk())

    grid  = case[:u].grid
    times = case[:u].times
    Nx, Ny, Nz = size(grid)

    VCCC = FieldTimeSeries{Center, Center, Center}(grid, times)
    VFCC = FieldTimeSeries{Face,   Center, Center}(grid, times)
    VCFC = FieldTimeSeries{Center, Face,   Center}(grid, times)
    VCCF = FieldTimeSeries{Center, Center, Face  }(grid, times)
    GC.gc()

    Nt = length(case[:u])

    params = Oceananigans.Utils.KernelParameters(0:Nx+1, 1:1, 0:Nz+1)
    _compute_volumes_kernel! = Oceananigans.Utils.configure_kernel(CPU(), grid, params, _compute_volumes!)[1]

    for t in 1:Nt
        _compute_volumes_kernel!(VCCC[t], VFCC[t], VCFC[t], VCCF[t], grid, case[:η][t])
    end
    
    case[:VCCC] = VCCC
    case[:VFCC] = VFCC
    case[:VCCF] = VCCF
    GC.gc()

    case[:Abx] = FieldTimeSeries(path, "Abx")
    case[:Abz] = FieldTimeSeries(path, "Abz")

    case[:abx] = [sum(case[:Abx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:abz] = [sum(case[:Abz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]
    case[:abt] = case[:abx] .+ case[:abz]

    case[:Gbx] = FieldTimeSeries(path, "Gbx")
    case[:Gbz] = FieldTimeSeries(path, "Gbz")
    GC.gc()

    case[:gbx] = [sum(case[:Gbx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:gbz] = [sum(case[:Gbz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]
    case[:gbt] = case[:gbx] .+ case[:gbz] 
    GC.gc()

    κb = FieldTimeSeries{Center, Center, Center}(grid, case[:Abx].times)

    for t in 1:Nt
        set!(κb[t], Diffusivity(case, t, :b))
    end

    case[:κb]  = κb
    case[:KE]  = [sum(u2(case, i))  + sum(w2(case, i))  for i in 1:Nt] ./ [sum(case[:VCCC][i]) for i in 1:Nt] ./ 2
    case[:MKE] = [sum(um2(case, i)) + sum(wm2(case, i)) for i in 1:Nt] ./ [sum(mean(case[:VCCC][i], dims=1)) for i in 1:Nt] ./ 2
    case[:η2]  = [mean(case[:η][i]^2) for i in 1:Nt]

    GC.gc()
    
    EDIAG = compute_rpe_density(case)

    case[:RPE] = EDIAG.rpe
    case[:APE] = EDIAG.ape

    return case
end

@inline function _diffusivity(i, j, k, grid, Ax, Az, Gz)
    ax = ℑxᶜᵃᵃ(i, j, k, grid, Ax)
    az = ℑzᵃᵃᶜ(i, j, k, grid, Az)
    gz = ℑzᵃᵃᶜ(i, j, k, grid, Gz)
    return - 0.5 * (ax + az) / gz
end

"""
    Diffusivity(case, t, var)

Compute the effective numerical diffusivity from variance dissipation diagnostics.

$(SIGNATURES)

# Arguments
- `case`: Case dictionary containing dissipation and gradient fields
- `t`: Time index
- `var`: Variable symbol (e.g., `:b` for buoyancy)

# Returns
- `KernelFunctionOperation` representing the numerical diffusivity field

The diffusivity is computed as `κ = -0.5 * (Ax + Az) / Gz`, where `Ax` and `Az` are the
advective dissipation rates and `Gz` is the vertical gradient squared. This provides a
local measure of numerical mixing, as described in the paper.
"""
function Diffusivity(case, t, var) 
    Ax = case[Symbol(:A, var, :x)][t]
    Az = case[Symbol(:A, var, :z)][t]
    Gz = case[Symbol(:G, var, :z)][t]
    grid = case[:u].grid
    return KernelFunctionOperation{Center, Center, Center}(_diffusivity, grid, Ax, Az, Gz)
end

