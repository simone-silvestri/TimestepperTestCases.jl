using Oceananigans.AbstractOperations: grid_metric_operation
using Oceananigans.Utils: launch!
using Oceananigans.Grids: architecture, znode
using Oceananigans.Architectures: device, on_architecture
using Oceananigans.Fields: default_indices
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index, @unroll

MetricField(loc, grid, metric; indices = default_indices(3)) = compute!(Field(grid_metric_operation(loc, metric, grid); indices))
@inline _density_operation(i, j, k, grid, b, ρ₀, g) = ρ₀ * (1 - b[i, j, k] / g)

DensityOperation(b; ρ₀ = 1000.0, g = 9.80655) = 
    KernelFunctionOperation{Center, Center, Center}(_density_operation, b.grid, b, ρ₀, g)

DensityField(b::Field; ρ₀ = 1000.0, g = 9.80655) = compute!(Field(DensityOperation(b; ρₒ, g)))

function compute_rpe_density(b::Field)
    ze = calculate_z★_diagnostics(b)
    εe = CenterField(b.grid) 
    αe = CenterField(b.grid) 
    zh = HeightField(ze.grid)

    ρ = DensityOperation(b)
    set!(εe, ze * ρ)
    set!(αe, (zh - ze) * ρ)

    return (; ze, εe, αe)
end

function HeightField(grid, loc=(Center(), Center(), Center()))  

    zf = Field(loc, grid)
    Lz = grid.Lz

    for k in 1:size(zf, 3)
        interior(zf, :, :, k) .= Lz + znode(k, grid, loc[3])
    end

    return zf
end

function calculate_z★_diagnostics(b::FieldTimeSeries; path = nothing)

    times = b.times

    if path isa Nothing
        path = b.path
    end

    vol = VolumeField(b.grid)
    z★  = FieldTimeSeries{Center, Center, Center}(b.grid, b.times; backend = OnDisk(), path, name = "z★")

    total_area = sum(AreaField(b.grid))
    
    z★t = CenterField(b.grid)

    for iter in 1:length(times)
        @info "time $iter of $(length(times))"

        calculate_z★!(z★t, b[iter], vol, total_area)
        set!(z★, z★t, iter)
    end
        
    return z★
end

function calculate_z★_diagnostics(b::Field)

    vol = VolumeField(b.grid)
    z★  = CenterField(b.grid)
    total_area = sum(AreaField(b.grid))
    calculate_z★!(z★, b, vol, total_area)
        
    return z★
end

function calculate_z★_diagnostics(b::FieldTimeSeries, i)

    times = b.times

    vol = VolumeField(b.grid)
    z★  = similar(b[1])

    total_area = sum(AreaField(b.grid))

    @info "time $i of $(length(times))"
    calculate_z★!(z★, b[i], vol, total_area)
        
    return z★
end

function calculate_z★!(z★::Field, b::Field, vol, total_area)
    grid = b.grid
    arch = architecture(grid)

    b_arr = Array(interior(b))[:]
    v_arr = Array(interior(vol))[:]

    perm           = sortperm(b_arr)
    sorted_b_field = b_arr[perm]
    sorted_v_field = v_arr[perm]
    integrated_v   = cumsum(sorted_v_field)    

    sorted_b_field = on_architecture(arch, sorted_b_field)
    integrated_v   = on_architecture(arch, integrated_v)

    launch!(arch, grid, :xyz, _calculate_z★, z★, b, sorted_b_field, integrated_v)
    
    z★ ./= total_area

    return nothing
end

@kernel function _calculate_z★(z★, b, b_sorted, integrated_v)
    i, j, k = @index(Global, NTuple)
    bl  = b[i, j, k]
    i₁  = searchsortedfirst(b_sorted, bl)
    z★[i, j, k] = integrated_v[i₁] 
end

function calculate_Γ²_diagnostics(z★::FieldTimeSeries, b::FieldTimeSeries; ρ₀ = 1000.0, g = 9.80655)
    
    times = b.times

    Γ²  = FieldTimeSeries{Center, Center, Center}(b.grid, b.times)

    for iter in 1:length(times)
        @info "time $iter of $(length(times))"

        ρ = DensityField(b[iter]; ρ₀, g)

        calculate_Γ²!(Γ²[iter], z★[iter], ρ)
    end
         
    return Γ²
end

function calculate_Γ²!(Γ², z★, ρ)
    grid = ρ.grid
    arch = architecture(grid)

    perm   = sortperm(Array(interior(z★))[:])

    ρ_arr  = (Array(interior(ρ))[:])[perm]
    z★_arr = (Array(interior(z★))[:])[perm]

    ρ_arr  = on_architecture(architecture(grid), ρ_arr)
    z★_arr = on_architecture(architecture(grid), z★_arr)
    launch!(arch, grid, :xyz, _calculate_Γ², Γ², z★, z★_arr, ρ_arr, grid)

    return nothing
end

@kernel function _calculate_Γ²(Γ², z★, z★_arr, ρ_arr, grid)
    i, j, k = @index(Global, NTuple)

    Nint = 10.0
     
    Γ²[i, j, k] = 0.0
         
    z_local  = znode(Center(), k, grid) + grid.Lz
    z★_local = z★[i, j, k] 
    Δz       = - (z_local - z★_local) / Nint
    zrange   = z_local:Δz:z★_local

    for z in zrange
        Γ²[i, j, k] += Δz * linear_interpolate(z★_arr, ρ_arr, z)
    end
end

