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
    case[:KE]  = [sum(u2(case, i))  + sum(w2(case, i))  for i in 1:Nt] ./ [sum(case[:VCCC][i]) for i in 1:Nt]
    case[:MKE] = [sum(um2(case, i)) + sum(wm2(case, i)) for i in 1:Nt] ./ [sum(mean(case[:VCCC][i], dims=1)) for i in 1:Nt]
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

function Diffusivity(case, t, var) 
    Ax = case[Symbol(:A, var, :x)][t]
    Az = case[Symbol(:A, var, :z)][t]
    Gz = case[Symbol(:G, var, :z)][t]
    grid = case[:u].grid
    return KernelFunctionOperation{Center, Center, Center}(_diffusivity, grid, Ax, Az, Gz)
end

