# # Internal tide: numerical mixing
#
# This example reads the internal-tide outputs and computes the horizontally averaged numerical diffusivity,
# its spatial and temporal averages, and differences relative to WRK3-SE. The stored dissipation and squared
# gradient fields include cell volumes; the existing loader performs the staggered-grid interpolation.
# Pass a local directory as the first command-line argument. Without it, the repository's internal_tide/
# directory is checked before the versioned download location in datasets.toml.

using TimestepperTestCases
using Oceananigans
using Statistics
include("data_sources.jl") #hide
include("compact_data.jl") #hide
include("plotting.jl") #hide

directory = isempty(ARGS) ? nothing : first(ARGS)
variants = discretizations()
profiles = NamedTuple[]

# Each case is reduced before the next one is loaded. The diffusivity is first computed locally as
# κ = −(Pₓ + P_z)/(2|∂z b|²), then averaged horizontally. The time average uses trapezoidal integration,
# so unequal output intervals are weighted by their duration.

for variant in variants
    filename = "internal_tide_$(variant.timestepper)_$(variant.label).jld2"
    data = compact_dataset("internal_tide", filename, prepare_internal_tide; directory)
    t, z = data.t, data.z
    P = data.fields["Abx"] + data.fields["Abz"]
    G = data.fields["Gbz"]
    κ = zeros(length(z), length(t))
    for n in eachindex(t), k in eachindex(z)
        wet = data.wet[:, k]
        κ[k, n] = mean(-P[wet, k, n] ./ (2G[wet, k, n]))
    end
    κᵗ = vec(sum((κ[:, 1:end-1] + κ[:, 2:end]) .* reshape(diff(t), 1, :) / 2; dims=2)) / (t[end] - t[1])
    κᶻ = vec(sum(κ .* data.Δz; dims=1)) / sum(data.Δz)
    push!(profiles, (; label=variant.label, t, z, κ, κᵗ, κᶻ))
    data = nothing
    GC.gc()
end

reference = only(filter(profile -> profile.label == "WRK3-SE", profiles))
differences = [(; profile..., Δκ=profile.κ - reference.κ) for profile in profiles if profile.label in ("AB2-SE", "MRK4-SE", "WRK3-UP")]

plot_internal_tide(profiles, differences) #hide
# ![Internal-tide numerical diffusivity](internal_tide_diffusivity.png)
