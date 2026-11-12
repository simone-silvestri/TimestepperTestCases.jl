using Oceananigans.Operators

"""
    load_channel(folder, case_number; arch=CPU())

Load channel snapshots and return velocity, buoyancy, free surface, cell volumes, and energy time series.
Saved interior values are used directly for the diagnostics.
"""
function load_channel(folder, case_number; arch = CPU())
    path = joinpath(folder, "snapshots_$(case_number).jld2")
    case = Dict()

    case[:u] = FieldTimeSeries(path, "u"; architecture=arch, backend=OnDisk(), boundary_conditions=nothing)
    case[:v] = FieldTimeSeries(path, "v"; architecture=arch, backend=OnDisk(), boundary_conditions=nothing)
    case[:w] = FieldTimeSeries(path, "w"; architecture=arch, backend=OnDisk(), boundary_conditions=nothing)
    case[:b] = FieldTimeSeries(path, "b"; architecture=arch, backend=OnDisk(), boundary_conditions=nothing)
    case[:η] = FieldTimeSeries(path, "η"; architecture=arch, backend=OnDisk(), boundary_conditions=nothing)

    grid = case[:u].grid
    Nx, Ny, Nz = size(grid)

    VCCC = FieldTimeSeries{Center, Center, Center}(grid, case[:u].times)
    VFCC = FieldTimeSeries{Face,   Center, Center}(grid, case[:u].times)
    VCFC = FieldTimeSeries{Center, Face,   Center}(grid, case[:u].times)
    VCCF = FieldTimeSeries{Center, Center, Face  }(grid, case[:u].times)
    GC.gc()

    Nt = length(case[:u])

    params = Oceananigans.Utils.KernelParameters(0:Nx+1, 0:Ny+1, 0:Nz+1)
    _compute_volumes_kernel! = Oceananigans.Utils.configure_kernel(arch, grid, params, _compute_volumes!)[1]

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
    case[:KE]  = [KineticEnergy(case, t) for t in 1:Nt]
    case[:MKE] = [MeanKineticEnergy(case, t) for t in 1:Nt]
    
    case[:η2]  = [mean(case[:η][i]^2) for i in 1:Nt]
    GC.gc()

    @info "Computing Potential Energy"
    EDIAG = compute_rpe_density_two(case)

    case[:RPE] = EDIAG.rpe
    case[:APE] = EDIAG.ape

    return case
end

@inline volume_weighted_square(i, j, k, grid, φ, V) = @inbounds φ[i, j, k]^2 * V[i, j, k]

"""
    KineticEnergy(case, i)

Return the volume-averaged kinetic energy, `(u² + v² + w²) / 2`, in m² s⁻² at snapshot `i`.
"""
function KineticEnergy(case, i)
    Vccc = case[:VCCC][i]
    Vfcc = case[:VFCC][i]
    Vcfc = case[:VCFC][i]
    Vccf = case[:VCCF][i]
    u    = case[:u][i]
    v    = case[:v][i]
    w    = case[:w][i]
    grid = u.grid

    u2 = KernelFunctionOperation{Face, Center, Center}(volume_weighted_square, grid, u, Vfcc)
    v2 = KernelFunctionOperation{Center, Face, Center}(volume_weighted_square, grid, v, Vcfc)
    w2 = KernelFunctionOperation{Center, Center, Face}(volume_weighted_square, grid, w, Vccf)
    
    return (sum(u2) + sum(v2) + sum(w2)) / (2sum(Vccc))
end

"""
    MeanKineticEnergy(case, i)

Return the volume-averaged kinetic energy of the zonal-mean velocity in m² s⁻² at snapshot `i`.
"""
function MeanKineticEnergy(case, i)
    Vccc = mean(case[:VCCC][i], dims=1)
    Vfcc = mean(case[:VFCC][i], dims=1)
    Vcfc = mean(case[:VCFC][i], dims=1)
    Vccf = mean(case[:VCCF][i], dims=1)
    u    = mean(case[:u][i], dims=1)
    v    = mean(case[:v][i], dims=1)
    w    = mean(case[:w][i], dims=1)
    grid = u.grid

    u2 = KernelFunctionOperation{Nothing, Center, Center}(volume_weighted_square, grid, u, Vfcc)
    v2 = KernelFunctionOperation{Nothing, Face,   Center}(volume_weighted_square, grid, v, Vcfc)
    w2 = KernelFunctionOperation{Nothing, Center, Face}(  volume_weighted_square, grid, w, Vccf)
    
    return (sum(u2) + sum(v2) + sum(w2)) / (2sum(Vccc))
end
