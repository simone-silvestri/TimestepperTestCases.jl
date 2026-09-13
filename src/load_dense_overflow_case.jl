using Oceananigans.ImmersedBoundaries: immersed_cell

"""
    load_dense_overflow(folder, free_surface_name, timestepper)

Open a dense overflow run and reduce it to the reference potential energy the case exists to measure.

$(SIGNATURES)

# Returns
- Dictionary with the fields `:u, :w, :T, :S, :b, :η`, the cell volumes `:VCCC`, the buoyancy variance
  dissipation `:Abx, :Abz` and its squared gradients `:Gbx, :Gbz`, the times `:times`, and the volume-averaged
  `:rpe` and `:ape` time series

The tracers and `b` are written volume weighted, so they are divided by `:VCCC` before use. The volumes are
rebuilt from `η` at every snapshot, the vertical coordinate following the free surface.
"""
function load_dense_overflow(folder, free_surface_name, timestepper)
    path = joinpath(folder, "dense_overflow_" * free_surface_name * "_" * timestepper * ".jld2")
    case = Dict()

    for name in (:u, :w, :T, :S, :b, :Gbx, :Gbz, :Abx, :Abz)
        case[name] = FieldTimeSeries(path, string(name); backend = OnDisk())
    end

    case[:η] = FieldTimeSeries(path, "η")
    fill_halo_regions!(case[:η])

    grid  = case[:u].grid
    times = case[:u].times
    case[:times] = times

    VCCC = FieldTimeSeries{Center, Center, Center}(grid, times)
    VFCC = FieldTimeSeries{Face,   Center, Center}(grid, times)
    VCFC = FieldTimeSeries{Center, Face,   Center}(grid, times)
    VCCF = FieldTimeSeries{Center, Center, Face  }(grid, times)

    for t in eachindex(times)
        launch!(CPU(), grid, :xyz, _compute_volumes!, VCCC[t], VFCC[t], VCFC[t], VCCF[t], grid, case[:η][t])
    end

    fill_halo_regions!(VCCC)
    case[:VCCC] = VCCC

    case[:rpe], case[:ape] = dense_overflow_rpe(case)

    return case
end

"""
    dense_overflow_rpe(case)

Volume-averaged reference and available potential energy of every snapshot, as `(rpe, ape)` [m² s⁻²].

$(SIGNATURES)

The dry cells are written as zeros, so the re-sorting is given the wet volume alone: a dry cell enters the
sort at `b = 0` carrying no volume, and contributes to neither the cumulative distribution nor the integral.

The weights are the volumes of the first snapshot, the bathymetry alone. The basin is not prismatic, so
`dz★/dV` varies by a factor of five over the column and the centimetres of free-surface breathing would move
`RPE` further than the mixing being measured does.
"""
function dense_overflow_rpe(case)
    grid = case[:b].grid
    Nx, Ny, Nz = size(grid)
    wet = [!immersed_cell(i, j, k, grid) for i in 1:Nx, j in 1:Ny, k in 1:Nz]

    b = Field{Center, Center, Center}(grid)
    V = Field{Center, Center, Center}(grid)

    interior(V) .= interior(case[:VCCC][1]) .* wet

    energies = map(eachindex(case[:times])) do t
        interior(b) .= interior(case[:b][t]) ./ interior(case[:VCCC][t])

        energy = compute_rpe_density(b, V)
        (sum(energy.εe * V) / sum(V), sum(energy.αe * V) / sum(V))
    end

    return first.(energies), last.(energies)
end

"""
    dense_overflow_diffusivity(case)

Effective diapycnal diffusivity [m² s⁻¹] of the run, `κ = (dRPE/dt) / N²`, with `N²` the equivalent
stratification of `dense_overflow_stability_parameters`.

$(SIGNATURES)

The configuration carries no tracer diffusivity, so the whole rise of the reference potential energy is
spurious and `κ` measures the discretization. The rate is the least-squares slope of `rpe` against `times`.
"""
function dense_overflow_diffusivity(case)
    t, rpe = case[:times], case[:rpe]
    t̄, r̄ = mean(t), mean(rpe)
    rate = sum((t .- t̄) .* (rpe .- r̄)) / sum((t .- t̄).^2)

    return rate / dense_overflow_stability_parameters().N²
end
