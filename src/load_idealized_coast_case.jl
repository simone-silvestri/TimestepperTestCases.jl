u2(case, i) = (case[:u][i])^2 / case[:VFCC][i]
v2(case, i) = (case[:v][i])^2 / case[:VCFC][i]
w2(case, i) = (case[:w][i])^2 / case[:VCCF][i]

um2(case, i) = (mean(case[:u][i], dims=1))^2 / mean(case[:VFCC][i], dims=1)
vm2(case, i) = (mean(case[:v][i], dims=1))^2 / mean(case[:VCFC][i], dims=1)
wm2(case, i) = (mean(case[:w][i], dims=1))^2 / mean(case[:VCCF][i], dims=1)

Sm(case, i)   = mean(case[:S][i] / case[:VCCC][i], dims=1)

function load_idealized_coast(folder, closure, suffix, timestepper)
    path = folder * "idealized_coast_" * suffix * timestepper * "_" * closure * ".jld2"
    case = Dict()

    case[:u] = FieldTimeSeries(path, "u") 
    case[:v] = FieldTimeSeries(path, "v")
    case[:w] = FieldTimeSeries(path, "w")
    case[:T] = FieldTimeSeries(path, "T")
    case[:S] = FieldTimeSeries(path, "S")
    case[:b] = FieldTimeSeries(path, "b")
    case[:η] = FieldTimeSeries(path, "η")

    grid = case[:u].grid
    Nx, Ny, Nz = size(grid)

    VCCC = FieldTimeSeries{Center, Center, Center}(grid, case[:u].times)
    VFCC = FieldTimeSeries{Face,   Center, Center}(grid, case[:u].times)
    VCFC = FieldTimeSeries{Center, Face,   Center}(grid, case[:u].times)
    VCCF = FieldTimeSeries{Center, Center, Face  }(grid, case[:u].times)

    Nt = length(case[:u])

    params = Oceananigans.Utils.KernelParameters(0:Nx+1, 0:Ny+1, 0:Nz+1)
    _compute_volumes_kernel! = Oceananigans.Utils.configure_kernel(CPU(), grid, params, _compute_volumes!)[1]

    for t in 1:Nt
        @info "Computing volumes for time step $t / $Nt"
        _compute_volumes_kernel!(VCCC[t], VFCC[t], VCFC[t], VCCF[t], grid, case[:η][t])
    end
    
    case[:VCCC] = VCCC
    case[:VFCC] = VFCC
    case[:VCFC] = VCFC
    case[:VCCF] = VCCF
    GC.gc()

    case[:ASx] = FieldTimeSeries(path, "ASx")
    case[:ASy] = FieldTimeSeries(path, "ASy")
    case[:ASz] = FieldTimeSeries(path, "ASz")

    case[:aSx] = [sum(case[:ASx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:aSy] = [sum(case[:ASy][i]) for i in 1:Nt] ./ [sum(case[:VCFC][i]) for i in 1:Nt]
    case[:aSz] = [sum(case[:ASz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]
    case[:aSt] = case[:aSx] .+ case[:aSy] .+ case[:aSz]

    case[:GSx] = FieldTimeSeries(path, "GSx")
    case[:GSy] = FieldTimeSeries(path, "GSy")
    case[:GSz] = FieldTimeSeries(path, "GSz")
    case[:DSz] = FieldTimeSeries(path, "DSz")

    case[:Gbx] = FieldTimeSeries(path, "Gbx")
    case[:Gby] = FieldTimeSeries(path, "Gby")
    case[:Gbz] = FieldTimeSeries(path, "Gbz")
    
    case[:gSx] = [sum(case[:GSx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:gSy] = [sum(case[:GSy][i]) for i in 1:Nt] ./ [sum(case[:VCFC][i]) for i in 1:Nt]
    case[:gSz] = [sum(case[:GSz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]
    case[:dSz] = [sum(case[:DSz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]

    case[:gbx] = [sum(case[:Gbx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:gby] = [sum(case[:Gby][i]) for i in 1:Nt] ./ [sum(case[:VCFC][i]) for i in 1:Nt]
    case[:gbz] = [sum(case[:Gbz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]
    
    case[:gSt] = case[:gSx] .+ case[:gSy] .+ case[:gSz]
    case[:gbt] = case[:gbx] .+ case[:gby] .+ case[:gbz]
    GC.gc()

    # case[:KE]  = [sum(u2(case, i))  + sum(v2(case, i))  + sum(w2(case, i))  for i in 1:Nt] ./ [sum(case[:VCCC][i]) for i in 1:Nt]
    # case[:MKE] = [sum(um2(case, i)) + sum(vm2(case, i)) + sum(wm2(case, i)) for i in 1:Nt] ./ [sum(mean(case[:VCCC][i], dims=1)) for i in 1:Nt]
    # case[:η2]  = [mean(case[:η][i]^2) for i in 1:Nt]

    GC.gc()
    EDIAG = compute_rpe_density(case)
    case[:RPE] = EDIAG.rpe
    case[:APE] = EDIAG.ape

    # κS = FieldTimeSeries{Center, Center, Center}(grid, case[:ASx].times)
    # κT = FieldTimeSeries{Center, Center, Center}(grid, case[:ASx].times)

    # for t in 1:Nt
    #     @info "Computing diffusivities for time step $t / $Nt"
    #     set!(κS[t], @at((Center, Center, Center),  (case[:ASx][t] + case[:ASy][t] + case[:ASz][t]) / 2 / (case[:GSx][t] + case[:GSy][t] + case[:GSz][t])))
    #     set!(κT[t], @at((Center, Center, Center),  (case[:ATx][t] + case[:ATy][t] + case[:ATz][t]) / 2 / (case[:GTx][t] + case[:GTy][t] + case[:GTz][t])))
    # end

    # case[:κS] = κS
    # case[:κT] = κT

    return case
end
