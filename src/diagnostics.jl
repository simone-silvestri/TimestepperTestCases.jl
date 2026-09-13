using Oceananigans.Utils: launch!
using Oceananigans.Grids: architecture, znode, znodes
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

@inline _density_operation(i, j, k, grid, b, ρ₀, g) = ρ₀ * (1 - b[i, j, k] / g)

DensityOperation(b; ρ₀ = 1000.0, g = 9.80655) = KernelFunctionOperation{Center, Center, Center}(_density_operation, b.grid, b, ρ₀, g)
DensityField(b::Field; ρ₀ = 1000.0, g = 9.80655) = compute!(Field(DensityOperation(b; ρₒ, g)))

"""
    compute_rpe_density(case::Dict)

Compute reference potential energy (RPE) and available potential energy (APE) time series from a case dictionary.

$(SIGNATURES)

# Arguments
- `case`: Dictionary containing `:b` (buoyancy FieldTimeSeries) and `:VCCC` (volume FieldTimeSeries)

# Returns
- Named tuple with `rpe` and `ape` arrays containing volume-averaged RPE and APE at each time step

RPE represents the minimum potential energy achievable by adiabatic re-sorting of the density field,
while APE is the difference between actual and reference potential energy. These diagnostics are
used to quantify irreversible mixing, as described in the paper.
"""
function compute_rpe_density(case::Dict)

    rpe = []
    ape = []

    for t in 1:length(case[:b])
        @info "time $t of $(length(case[:b]))"
        push!(rpe, sum(compute_rpe_density(case[:b][t] / case[:VCCC][t], case[:VCCC][t]).εe * on_architecture(CPU(), case[:VCCC][t])) / sum(case[:VCCC][t]))
        push!(ape, sum(compute_rpe_density(case[:b][t] / case[:VCCC][t], case[:VCCC][t]).αe * on_architecture(CPU(), case[:VCCC][t])) / sum(case[:VCCC][t]))
    end

    return (; rpe, ape)
end

"""
    compute_rpe_density_two(case::Dict)

Compute RPE and APE time series using buoyancy directly (without volume weighting).

$(SIGNATURES)

# Arguments
- `case`: Dictionary containing `:b` (buoyancy FieldTimeSeries) and `:VCCC` (volume FieldTimeSeries)

# Returns
- Named tuple with `rpe` and `ape` arrays

This variant computes RPE/APE using the buoyancy field directly, suitable for cases where
buoyancy is already volume-weighted or when using different volume conventions.
"""
function compute_rpe_density_two(case::Dict)

    rpe = []
    ape = []

    for t in 1:length(case[:b])
        @info "time $t of $(length(case[:b]))"
        push!(rpe, sum(compute_rpe_density(case[:b][t], case[:VCCC][t]).εe * on_architecture(CPU(), case[:VCCC][t])) / sum(case[:VCCC][t]))
        push!(ape, sum(compute_rpe_density(case[:b][t], case[:VCCC][t]).αe * on_architecture(CPU(), case[:VCCC][t])) / sum(case[:VCCC][t]))
    end

    return (; rpe, ape)
end

compute_rpe_density(b::Oceananigans.AbstractOperations.AbstractOperation, vol) = compute_rpe_density(Field(b), vol)

"""
    compute_rpe_density(b::Field, vol)

Compute RPE and APE density fields from buoyancy and volume fields.

$(SIGNATURES)

# Arguments
- `b`: Buoyancy field
- `vol`: Volume field

# Returns
- Named tuple with `ze` (re-sorted height), `εe` (RPE density), and `αe` (APE density) fields

This function computes the re-sorted height `z★` by sorting buoyancy and computing the cumulative
volume distribution. RPE density is `ρ * z★` and APE density is `ρ * (z - z★)`, where `ρ` is
computed from buoyancy using the linear equation of state.
"""
function compute_rpe_density(b::Field, vol)
    bcpu = on_architecture(CPU(), b)
    volcpu = on_architecture(CPU(), vol)
    ze = calculate_z★_diagnostics(bcpu, volcpu)
    εe = CenterField(bcpu.grid) 
    αe = CenterField(bcpu.grid) 
    zh = HeightField(ze.grid)

    ρ = DensityOperation(bcpu)
    set!(εe, ze * ρ)
    set!(αe, (zh - ze) * ρ)

    return (; ze, εe, αe)
end

"""
    calculate_z★_diagnostics(b::Field, vol)

Calculate the re-sorted height field `z★` from buoyancy and volume fields.

$(SIGNATURES)

# Arguments
- `b`: Buoyancy field
- `vol`: Volume field

# Returns
- `z★` field representing the height each fluid parcel would occupy after adiabatic re-sorting

The re-sorted height is computed by sorting the buoyancy field and computing cumulative volume
distribution, providing a reference state for potential energy calculations.
"""
function calculate_z★_diagnostics(b::Field, vol)

    z★ = CenterField(b.grid)
    calculate_z★!(z★, b, vol)
        
    return z★
end

function calculate_z★_diagnostics(b::FieldTimeSeries, i)

    times = b.times

    vol = VolumeField(b.grid)
    z★  = similar(b[1])

    @info "time $i of $(length(times))"
    calculate_z★!(z★, b[i], vol)
        
    return z★
end

# Height at which the volume below equals `Vc`, by inverting the cumulative volume `volume_below` of the faces
# `zf`. A flat-bottomed column of the same total volume would place the sorted state above the real one
# wherever the basin is not prismatic, and the available energy would come out negative.
@inline function sorted_height(volume_below, zf, Vc)
    k  = clamp(searchsortedlast(volume_below, Vc), 1, length(zf) - 1)
    ΔV = volume_below[k+1] - volume_below[k]
    θ  = ΔV > 0 ? (Vc - volume_below[k]) / ΔV : zero(Vc)

    return zf[k] + θ * (zf[k+1] - zf[k])
end

function calculate_z★!(z★::Field, b::Field, vol)
    b_arr = Array(interior(b))[:]
    v_arr = Array(interior(vol))[:]

    valid = (b_arr .!= 0) .& (!).(isnan.(b_arr))

    if !any(valid)
        @warn "calculate_z★!: no valid buoyancy cells, returning NaN z★"
        fill!(z★, convert(eltype(z★), NaN))
        return nothing
    end

    wet_volume   = reshape(v_arr .* valid, size(interior(vol)))
    volume_below = [0.0; cumsum(vec(sum(wet_volume, dims = (1, 2))))]
    zf           = Array(znodes(b.grid, Face()))

    cells = findall(valid)
    perm  = sortperm(b_arr[cells])

    # each parcel takes the cumulative volume at its own rank in the sorted order: looking z★ up by buoyancy
    # gives every member of a block of equal buoyancy the volume of the whole block, which biases z★ high and
    # makes the reference potential energy fall as the ties break up
    sorted_volume = v_arr[cells][perm]
    midpoint = cumsum(sorted_volume) .- sorted_volume ./ 2

    sorted = zeros(eltype(z★), length(b_arr))
    sorted[cells[perm]] .= sorted_height.(Ref(volume_below), Ref(zf), midpoint)

    interior(z★) .= reshape(sorted, size(interior(z★)))

    return nothing
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

