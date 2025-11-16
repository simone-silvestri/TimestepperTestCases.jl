using Oceananigans.AbstractOperations: grid_metric_operation
using Oceananigans.Utils: launch!
using Oceananigans.Grids: architecture, znode
using Oceananigans.Architectures: device, on_architecture
using Oceananigans.Fields: default_indices
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index

@kernel function _compute_volumes!(VCCC, VFCC, VCFC, VCCF, grid, η)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        VCCC[i, j, k] = Δxᶜᶜᶜ(i, j, k, grid) * Δyᶜᶜᶜ(i, j, k, grid) * Δrᶜᶜᶜ(i, j, k, grid) * column_depthᶜᶜᵃ(i, j, grid.Nz+1, grid, η) / static_column_depthᶜᶜᵃ(i, j, grid)
        VFCC[i, j, k] = Δxᶠᶜᶜ(i, j, k, grid) * Δyᶠᶜᶜ(i, j, k, grid) * Δrᶠᶜᶜ(i, j, k, grid) * column_depthᶠᶜᵃ(i, j, grid.Nz+1, grid, η) / static_column_depthᶠᶜᵃ(i, j, grid)
        VCFC[i, j, k] = Δxᶜᶠᶜ(i, j, k, grid) * Δyᶜᶠᶜ(i, j, k, grid) * Δrᶜᶠᶜ(i, j, k, grid) * column_depthᶜᶠᵃ(i, j, grid.Nz+1, grid, η) / static_column_depthᶜᶠᵃ(i, j, grid)
        VCCF[i, j, k] = Δxᶜᶜᶠ(i, j, k, grid) * Δyᶜᶜᶠ(i, j, k, grid) * Δrᶜᶜᶠ(i, j, k, grid) * column_depthᶜᶜᵃ(i, j, grid.Nz+1, grid, η) / static_column_depthᶜᶜᵃ(i, j, grid)
    end
end

HeightField(grid) = Field(KernelFunctionOperation{Center, Center, Center}(Oceananigans.Grids.znode, grid, Center(), Center(), Center()))
AreaField(grid)   = Field(KernelFunctionOperation{Center, Center, Nothing}(Oceananigans.Operators.Azᶜᶜᶜ, grid))

@inline _density_operation(i, j, k, grid, b, ρ₀, g) = ρ₀ * (1 - b[i, j, k] / g)

DensityOperation(b; ρ₀ = 1000.0, g = 9.80655) = KernelFunctionOperation{Center, Center, Center}(_density_operation, b.grid, b, ρ₀, g)
DensityField(b::Field; ρ₀ = 1000.0, g = 9.80655) = compute!(Field(DensityOperation(b; ρₒ, g)))

function compute_rpe_density(case::Dict)

    rpe = []
    ape = []

    for t in 1:length(case[:b])
        @info "time $t of $(length(case[:b]))"
        push!(rpe, sum(compute_rpe_density(case[:b][t] / case[:VCCC][t], case[:VCCC][t]).εe * case[:VCCC][t]) / sum(case[:VCCC][t]))
        push!(ape, sum(compute_rpe_density(case[:b][t] / case[:VCCC][t], case[:VCCC][t]).αe * case[:VCCC][t]) / sum(case[:VCCC][t]))
    end

    return (; rpe, ape)
end

function compute_rpe_density_two(case::Dict)

    rpe = []
    ape = []

    for t in 1:length(case[:b])
        @info "time $t of $(length(case[:b]))"
        push!(rpe, sum(compute_rpe_density(case[:b][t], case[:VCCC][t]).εe * case[:VCCC][t]) / sum(case[:VCCC][t]))
        push!(ape, sum(compute_rpe_density(case[:b][t], case[:VCCC][t]).αe * case[:VCCC][t]) / sum(case[:VCCC][t]))
    end

    return (; rpe, ape)
end

compute_rpe_density(b::Oceananigans.AbstractOperations.AbstractOperation, vol) = compute_rpe_density(Field(b), vol)

function compute_rpe_density(b::Field, vol)
    ze = calculate_z★_diagnostics(b, vol)
    εe = CenterField(b.grid) 
    αe = CenterField(b.grid) 
    zh = HeightField(ze.grid)

    ρ = DensityOperation(b)
    set!(εe, ze * ρ)
    set!(αe, (zh - ze) * ρ)

    return (; ze, εe, αe)
end

function calculate_z★_diagnostics(b::Field, vol)

    total_area = sum(AreaField(b.grid))
    z★ = CenterField(b.grid)
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

    valid_indices = (b_arr .!= 0) .& (!).(isnan.(b_arr))
    b_arr = b_arr[valid_indices]
    v_arr = v_arr[valid_indices]

    perm           = sortperm(b_arr)
    sorted_b_field = b_arr[perm]
    sorted_v_field = v_arr[perm]
    integrated_v   = cumsum(sorted_v_field)    

    launch!(arch, grid, :xyz, _calculate_z★, z★, b, sorted_b_field, integrated_v)
    
    z★ ./= total_area

    return nothing
end

@kernel function _calculate_z★(z★, b, b_sorted, integrated_v)
    i, j, k = @index(Global, NTuple)
    bl  = b[i, j, k]
    i₁  = searchsortedlast(b_sorted, bl)
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

