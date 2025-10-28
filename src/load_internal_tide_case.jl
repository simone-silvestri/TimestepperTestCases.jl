function load_internal_tide(folder, suffix, timestepper)
    path = folder * "internal_tide_" * suffix * timestepper * ".jld2"
    case = Dict()

    case[:u] = FieldTimeSeries(path, "u")
    case[:w] = FieldTimeSeries(path, "w")
    case[:b] = FieldTimeSeries(path, "b")
    case[:η] = FieldTimeSeries(path, "η")

    grid = case[:u].grid
    Nx, Ny, Nz = size(grid)

    VCCC = FieldTimeSeries{Center, Center, Center}(grid, case[:u].times)
    VFCC = FieldTimeSeries{Face,   Center, Center}(grid, case[:u].times)
    VCFC = FieldTimeSeries{Center, Face,   Center}(grid, case[:u].times)
    VCCF = FieldTimeSeries{Center, Center, Face  }(grid, case[:u].times)

    Nt = length(case[:u])

    params = Oceananigans.Utils.KernelParameters(0:Nx+1, 1:1, 0:Nz+1)
    _compute_volumes_kernel! = Oceananigans.Utils.configure_kernel(CPU(), grid, params, _compute_volumes!)[1]

    for t in 1:Nt
        @info "Computing volumes for time step $t / $Nt"
        _compute_volumes_kernel!(VCCC[t], VFCC[t], VCFC[t], VCCF[t], grid, case[:η][t])
    end
    
    case[:VCCC] = VCCC
    case[:VFCC] = VFCC
    case[:VCCF] = VCCF
    GC.gc()

    case[:Abx] = FieldTimeSeries(path, "Abx")
    case[:Abz] = FieldTimeSeries(path, "Abz")
    case[:Acx] = FieldTimeSeries(path, "Acx")
    case[:Acz] = FieldTimeSeries(path, "Acz")


    case[:abx] = [sum(case[:Abx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:abz] = [sum(case[:Abz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]
    case[:abt] = case[:abx] .+ case[:abz]

    case[:acx] = [sum(case[:Acx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:acz] = [sum(case[:Acz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]
    case[:act] = case[:acx] .+ case[:acz]

    case[:Gbx] = FieldTimeSeries(path, "Gbx")
    case[:Gbz] = FieldTimeSeries(path, "Gbz")
    case[:Gcx] = FieldTimeSeries(path, "Gcx")
    case[:Gcz] = FieldTimeSeries(path, "Gcz")
    
    case[:gbx] = [sum(case[:Gbx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:gbz] = [sum(case[:Gbz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]
    case[:gcx] = [sum(case[:Gcx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:gcz] = [sum(case[:Gcz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]

    case[:gbt] = case[:gbx] .+ case[:gbz] 
    case[:gct] = case[:gcx] .+ case[:gcz] 
    GC.gc()

    case[:KE]  = [sum(u2(case, i))  + sum(w2(case, i))  for i in 1:Nt] ./ [sum(case[:VCCC][i]) for i in 1:Nt]
    case[:MKE] = [sum(um2(case, i)) + sum(wm2(case, i)) for i in 1:Nt] ./ [sum(mean(case[:VCCC][i], dims=1)) for i in 1:Nt]
    case[:η2]  = [mean(case[:η][i]^2) for i in 1:Nt]

    GC.gc()
    EDIAG = compute_rpe_density(case)

    case[:RPE] = EDIAG.rpe
    case[:APE] = EDIAG.ape

    return case
end
