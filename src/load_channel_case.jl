using Oceananigans.Operators

u2a(case, i) = (case[:u][i])^2 * case[:VFCC][i]
v2a(case, i) = (case[:v][i])^2 * case[:VCFC][i]
w2a(case, i) = (case[:w][i])^2 * case[:VCCF][i]

um2a(case, i) = (mean(case[:u][i], dims=1))^2 * mean(case[:VFCC][i], dims=1)
vm2a(case, i) = (mean(case[:v][i], dims=1))^2 * mean(case[:VCFC][i], dims=1)
wm2a(case, i) = (mean(case[:w][i], dims=1))^2 * mean(case[:VCCF][i], dims=1)

"""
    load_channel(folder, case_number; arch)

Load and process channel simulation output data.

$(SIGNATURES)

# Arguments
- `folder`: Directory containing the simulation output files
- `case_number`: Case identifier string
- `arch`: Architecture to load data on (default: `CPU()`)

# Returns
- Dictionary containing velocity, buoyancy, volume, and energy diagnostics

This function loads channel simulation snapshots and computes kinetic energy, mean kinetic
energy, and potential energy diagnostics for equilibrated channel flow analysis.
"""
function load_channel(folder, case_number; arch = CPU())
    path = folder * "snapshots_$(case_number).jld2"
    @show path
    case = Dict()

    case[:u] = FieldTimeSeries(path, "u"; architecture=arch, backend=OnDisk())
    case[:v] = FieldTimeSeries(path, "v"; architecture=arch, backend=OnDisk())
    case[:w] = FieldTimeSeries(path, "w"; architecture=arch, backend=OnDisk())
    case[:b] = FieldTimeSeries(path, "b"; architecture=arch, backend=OnDisk())
    case[:η] = FieldTimeSeries(path, "η"; architecture=arch, backend=OnDisk())

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
    case[:KE]  = [sum(KineticEnergy(case, 1)) / sum(case[:VCCC][1])]
    case[:MKE] = [sum(MeanKineticEnergy(case, 1)) / sum(mean(case[:VCCC][1], dims=1))]

    for t in 2:Nt
        @info "computing index $t"
        push!(case[:KE] , KineticEnergy(case, t))
        push!(case[:MKE], MeanKineticEnergy(case, t))
    end
    
    case[:η2]  = [mean(case[:η][i]^2) for i in 1:Nt]
    GC.gc()

    @info "Computing Potential Energy"
    EDIAG = compute_rpe_density_two(case)

    case[:RPE] = EDIAG.rpe
    case[:APE] = EDIAG.ape

    return case
end

@inline _Vφ²(i, j, k, grid, φ, V) = @inbounds φ[i, j, k]^2 * V[i, j, k]

"""
    KineticEnergy(case, i)

Compute the volume-averaged kinetic energy at time index `i`.

$(SIGNATURES)

# Arguments
- `case`: Case dictionary containing velocity and volume fields
- `i`: Time index

# Returns
- Volume-averaged kinetic energy [m²/s²]: `(u² + v² + w²) / V`

Kinetic energy is computed using volume-weighted velocity components at their respective
grid locations (u at Face-Center-Center, v at Center-Face-Center, w at Center-Center-Face).
"""
function KineticEnergy(case, i)
    Vccc = case[:VCCC][i]
    Vfcc = case[:VFCC][i]
    Vcfc = case[:VCFC][i]
    Vccf = case[:VCCF][i]
    u = case[:u][i]
    v = case[:v][i]
    w = case[:w][i]
    grid = u.grid

    u2 = KernelFunctionOperation{Face, Center, Center}(_Vφ², grid, u, Vfcc)
    v2 = KernelFunctionOperation{Center, Face, Center}(_Vφ², grid, v, Vcfc)
    w2 = KernelFunctionOperation{Center, Center, Face}(_Vφ², grid, w, Vccf)
    
    return (sum(u2) + sum(v2) + sum(w2)) / sum(Vccc)
end

"""
    MeanKineticEnergy(case, i)

Compute the volume-averaged mean kinetic energy (zonally averaged) at time index `i`.

$(SIGNATURES)

# Arguments
- `case`: Case dictionary containing velocity and volume fields
- `i`: Time index

# Returns
- Volume-averaged mean kinetic energy [m²/s²]

Mean kinetic energy is computed from zonally averaged velocity components, representing
the energy in the mean flow rather than eddy kinetic energy.
"""
function MeanKineticEnergy(case, i)
    Vccc = mean(case[:VCCC][i], dims=1)
    Vfcc = mean(case[:VFCC][i], dims=1)
    Vcfc = mean(case[:VCFC][i], dims=1)
    Vccf = mean(case[:VCCF][i], dims=1)
    u = mean(case[:u][i], dims=1) 
    v = mean(case[:v][i], dims=1) 
    w = mean(case[:w][i], dims=1) 
    grid = u.grid

    u2 = KernelFunctionOperation{Nothing, Center, Center}(_Vφ², grid, u, Vfcc)
    v2 = KernelFunctionOperation{Nothing, Face,   Center}(_Vφ², grid, v, Vcfc)
    w2 = KernelFunctionOperation{Nothing, Center, Face}(  _Vφ², grid, w, Vccf)
    
    return (sum(u2) + sum(v2) + sum(w2)) / sum(Vccc)
end
