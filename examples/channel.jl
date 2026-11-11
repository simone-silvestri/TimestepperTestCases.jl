# # Equilibrated channel: numerical diffusivity and dissipation
#
# The input files contain the averaged buoyancy-variance tendencies and squared gradients. This example uses
# the final saved averaging window, recorded below for each case. Only the averages files are required;
# three-dimensional snapshots and restart files are unnecessary for these diagnostics.

using TimestepperTestCases
using Oceananigans
using Statistics
include("data_sources.jl") #src
include("compact_data.jl") #src
include("plotting.jl") #src

directory = isempty(ARGS) ? nothing : first(ARGS)
labels = ("AB2-SE", "WRK3-SE-SM05", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM", "WRK3-UP")
profiles = NamedTuple[]
sections = NamedTuple[]

# The flux components live on different faces. Interpolation to cell centers precedes their sum. Numerical
# diffusivity is the ratio of the horizontally averaged variance tendency to twice the averaged squared gradient.

for label in labels
    data = compact_dataset("channel_simulation", "averages_$label.jld2", prepare_channel; directory)
    @info "Channel averaging window" case=label final_year=data.time / (360 * 86400)
    Pₕ = data.fields["Abx"] + data.fields["Aby"]
    Pz = data.fields["Abz"]
    Gₕ = data.fields["Gbx"] + data.fields["Gby"]
    Gz = data.fields["Gbz"]
    κₕ = vec(-sum(Pₕ; dims=1) ./ (2sum(Gₕ; dims=1)))
    κ = vec(-sum(Pₕ + Pz; dims=1) ./ (2sum(Gₕ + Gz; dims=1)))
    κᵖ = vec(-sum(data.fields["Dbz"]; dims=1) ./ (2sum(Gz; dims=1)))
    push!(profiles, (; label, data.z, κₕ, κ, ratio=κᵖ ./ κ, data.time))
    push!(sections, (; label, data.y, data.z, Pₕ, Pz))
end

# The section differences retain the sign of the variance tendency and use WRK3-SE as the reference.
reference = only(filter(section -> section.label == "WRK3-SE", sections))
differences = [(; section..., ΔPₕ=section.Pₕ - reference.Pₕ, ΔPz=section.Pz - reference.Pz) for section in sections if section.label in ("AB2-SE", "MRK4-SE", "WRK3-UP")]

plot_channel(profiles, reference, differences) #src
# ![Channel numerical diffusivity](channel_diffusivity.png)
# ![Channel variance-tendency differences](channel_dissipation.png)
