# # Equilibrated channel: numerical diffusivity and dissipation
#
# The input files contain the averaged buoyancy-variance tendencies and squared gradients. This example uses
# the final saved averaging window, recorded below for each case. Only the averages files are required;
# three-dimensional snapshots and restart files are unnecessary for these diagnostics.

using TimestepperTestCases
using Oceananigans
using Statistics
include("data_sources.jl") #hide
include("plotting.jl") #hide

directory = isempty(ARGS) ? nothing : first(ARGS)
labels = ("AB2-SE", "WRK3-SE", "MRK4-SE", "WRK3-UP")
profiles = NamedTuple[]
sections = NamedTuple[]

# The flux components live on different faces. Interpolation to cell centers precedes their sum. Numerical
# diffusivity is the ratio of the horizontally averaged variance tendency to twice the averaged squared gradient.

for label in labels
    path = dataset_file("channel_simulation", "averages_$label.jld2"; directory)
    fields = Dict(name => FieldTimeSeries(path, name; backend=OnDisk()) for name in ("Abx", "Aby", "Abz", "Gbx", "Gby", "Gbz", "Dbz"))
    grid = fields["Abx"].grid
    n = length(fields["Abx"].times)
    centered = Dict{String, Array{Float64, 3}}()
    for name in keys(fields)
        field = Field{Center, Center, Center}(grid)
        set!(field, fields[name][n])
        centered[name] = Array(interior(field))
    end
    Pₕ = centered["Abx"] + centered["Aby"]
    Pz = centered["Abz"]
    Gₕ = centered["Gbx"] + centered["Gby"]
    Gz = centered["Gbz"]
    κₕ = vec(-sum(Pₕ; dims=(1, 2)) ./ (2sum(Gₕ; dims=(1, 2))))
    κ = vec(-sum(Pₕ + Pz; dims=(1, 2)) ./ (2sum(Gₕ + Gz; dims=(1, 2))))
    κᵖ = vec(-sum(centered["Dbz"]; dims=(1, 2)) ./ (2sum(Gz; dims=(1, 2))))
    z = collect(znodes(grid, Center()))
    y = collect(ynodes(grid, Center()))
    push!(profiles, (; label, z, κₕ, κ, ratio=κᵖ ./ κ, time=fields["Abx"].times[n]))
    push!(sections, (; label, y, z, Pₕ=dropdims(mean(Pₕ; dims=1); dims=1), Pz=dropdims(mean(Pz; dims=1); dims=1)))
end

# The section differences retain the sign of the variance tendency and use WRK3-SE as the reference.
reference = only(filter(section -> section.label == "WRK3-SE", sections))
differences = [(; section..., ΔPₕ=section.Pₕ - reference.Pₕ, ΔPz=section.Pz - reference.Pz) for section in sections if section.label != reference.label]

plot_channel(profiles, reference, differences) #hide
# ![Channel numerical diffusivity](channel_diffusivity.png)
# ![Channel variance-tendency differences](channel_dissipation.png)
