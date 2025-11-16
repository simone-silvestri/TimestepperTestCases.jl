using Oceananigans.Operators
import Oceananigans.BoundaryConditions: fill_halo_event!

const NoBCs = Union{Nothing, Missing, Tuple{Vararg{Nothing}}}

@inline fill_halo_event!(c, ::Nothing, bcs, loc, grid, args...; kwargs...) = nothing
@inline fill_halo_event!(c, kernel!, ::NoBCs, loc, grid, args...; kwargs...) = nothing

u2a(case, i) = (case[:u][i])^2 * case[:VFCC][i]
v2a(case, i) = (case[:v][i])^2 * case[:VCFC][i]
w2a(case, i) = (case[:w][i])^2 * case[:VCCF][i]

um2a(case, i) = (mean(case[:u][i], dims=1))^2 * mean(case[:VFCC][i], dims=1)
vm2a(case, i) = (mean(case[:v][i], dims=1))^2 * mean(case[:VCFC][i], dims=1)
wm2a(case, i) = (mean(case[:w][i], dims=1))^2 * mean(case[:VCCF][i], dims=1)

Sm(case, i)   = mean(case[:S][i] / case[:VCCC][i], dims=1)

function load_channel(folder, case_number)
    path = folder * "snapshots_$(case_number).jld2"
    @show path
    case = Dict()

    case[:u] = FieldTimeSeries(path, "u"; backend=OnDisk())
    case[:v] = FieldTimeSeries(path, "v"; backend=OnDisk())
    case[:w] = FieldTimeSeries(path, "w"; backend=OnDisk())
    case[:b] = FieldTimeSeries(path, "b"; backend=OnDisk())
    case[:η] = FieldTimeSeries(path, "η"; backend=OnDisk())

    grid = case[:u].grid
    Nx, Ny, Nz = size(grid)

    VCCC = FieldTimeSeries{Center, Center, Center}(grid, case[:u].times)
    VFCC = FieldTimeSeries{Face,   Center, Center}(grid, case[:u].times)
    VCFC = FieldTimeSeries{Center, Face,   Center}(grid, case[:u].times)
    VCCF = FieldTimeSeries{Center, Center, Face  }(grid, case[:u].times)
    GC.gc()

    Nt = length(case[:u])

    params = Oceananigans.Utils.KernelParameters(0:Nx+1, 0:Ny+1, 0:Nz+1)
    _compute_volumes_kernel! = Oceananigans.Utils.configure_kernel(CPU(), grid, params, _compute_volumes!)[1]

    for t in 1:Nt
        @info "Computing volumes $t of $Nt" 
        _compute_volumes_kernel!(VCCC[t], VFCC[t], VCFC[t], VCCF[t], grid, case[:η][t])
    end
    
    case[:VCCC] = VCCC
    case[:VFCC] = VFCC
    case[:VCFC] = VCFC
    case[:VCCF] = VCCF
    GC.gc()

    @info "Computing Kinetic Energy"
    case[:KE]  = [sum(KineticEnergy(case, 1)) / sum(case[:VCCC][1])]
    case[:MKE] = [sum(MeanKineticEnergy(case, 1)) / sum(mean(case[:VCCC][1], dims=1))]

    for t in 2:Nt
        @info "computing index $t"
        push!(case[:KE] , sum(KineticEnergy(case, t)) / sum(case[:VCCC][t]))
        push!(case[:MKE], sum(MeanKineticEnergy(case, t)) / sum(mean(case[:VCCC][t], dims=1)))
    end
    
    case[:η2]  = [mean(case[:η][i]^2) for i in 1:Nt]
    GC.gc()

    @info "Computing Potential Energy"
    EDIAG = compute_rpe_density_two(case)

    case[:RPE] = EDIAG.rpe
    case[:APE] = EDIAG.ape

    return case
end

@inline Vφ²(i, j, k, grid, φ, V) = @inbounds φ[i, j, k]^2 * V[i, j, k]

@inline function _kinetic_energy(i, j, k, grid, u, v, w, Vfcc, Vcfc, Vccf)
    u2n = ℑxᶜᵃᵃ(i, j, k, grid, Vφ², u, Vfcc)
    v2n = ℑyᵃᶜᵃ(i, j, k, grid, Vφ², v, Vcfc)
    w2n = ℑzᵃᵃᶜ(i, j, k, grid, Vφ², w, Vccf)
    return u2n + v2n + w2n
end

function KineticEnergy(case, i)
    Vccc = case[:VCCC][i]
    Vfcc = case[:VFCC][i]
    Vcfc = case[:VCFC][i]
    Vccf = case[:VCCF][i]
    u = case[:u][i]
    v = case[:v][i]
    w = case[:w][i]
    grid = u.grid

    return KernelFunctionOperation{Center, Center, Center}(_kinetic_energy, grid, u, v, w, Vfcc, Vcfc, Vccf)
end

function MeanKineticEnergy(case, i)
    Vccc = mean(case[:VCCC][i], dims=1)
    Vfcc = mean(case[:VFCC][i], dims=1)
    Vcfc = mean(case[:VCFC][i], dims=1)
    Vccf = mean(case[:VCCF][i], dims=1)
    u = mean(case[:u][i], dims=1) 
    v = mean(case[:v][i], dims=1) 
    w = mean(case[:w][i], dims=1) 
    grid = u.grid

    return KernelFunctionOperation{Nothing, Center, Center}(_kinetic_energy, grid, u, v, w, Vfcc, Vcfc, Vccf)
end
