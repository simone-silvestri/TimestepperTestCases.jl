# # Internal tide: numerical mixing
#
# This example reads the internal-tide outputs and computes the horizontally averaged numerical diffusivity,
# its spatial and temporal averages, and differences relative to WRK3-SE. The stored dissipation and squared
# gradient fields include cell volumes and are interpolated to cell centers before reduction.
# Pass a local directory as the first command-line argument. Without it, the repository's internal_tide/
# directory is checked before the versioned download location in datasets.toml.

using TimestepperTestCases
using Oceananigans
using Statistics
include("../shared/data_sources.jl") #src
include("../shared/compact_data.jl") #src
include("../shared/plotting.jl") #src
function plot_internal_tide(profiles, differences)
    figure       = Figure(size=(1100, 720))
    history      = Axis(figure[1, 1:2]; xlabel=L"t\;[\mathrm{days}]", ylabel=L"\langle\kappa\rangle\;[\mathrm{m^2\,s^{-1}}]")
    profile_axis = Axis(figure[1, 3]; xlabel=L"\overline{\kappa}\;[\mathrm{m^2\,s^{-1}}]", ylabel=L"z\;[\mathrm{km}]")
    for profile in profiles
        style = plot_style(profile.label)
        lines!(history, profile.t / 86400, profile.κᶻ; style..., label=profile.label)
        lines!(profile_axis, profile.κᵗ, profile.z / 1000; style...)
    end
    Legend(figure[2, 1:3], history; nbanks=4, framevisible=true)
    rowsize!(figure.layout, 2, Auto(false))
    rowgap!(figure.layout, 16)
    limit = 1.5e-5
    for (column, profile) in enumerate(differences)
        axis = Axis(figure[3, column]; title=profile.label * " − WRK3-SE", xlabel=L"t\;[\mathrm{days}]", ylabel=L"z\;[\mathrm{km}]")
        heatmap!(axis, profile.t / 86400, profile.z / 1000, transpose(profile.Δκ);
                 colormap=:balance, colorrange=(-limit, limit), colorscale=Makie.Symlog10(limit / 100))
    end
    Colorbar(figure[3, 4]; colormap=:balance, limits=(-limit, limit), scale=Makie.Symlog10(limit / 100), label=L"\Delta\kappa\;[\mathrm{m^2\,s^{-1}}]")
    return save_example(figure, "internal_tide_diffusivity")
end

directory = isempty(ARGS) ? nothing : first(ARGS)
variants = discretizations()
profiles = NamedTuple[]

# Numerical diffusivity is the absolute ratio of horizontal averages, with ε = 10⁻⁷ added to the
# volume-integrated squared-gradient denominator. A box average spans four depth indices and six time indices
# on each side. The time mean uses trapezoidal integration over the saved output times.

ϵκ = 1e-7

for variant in variants
    filename = "internal_tide_$(variant.timestepper)_$(variant.label).jld2"
    data     = compact_dataset("internal_tide", filename, prepare_internal_tide; directory)
    t, z = data.t, data.z
    P  = data.fields["Abx"] + data.fields["Abz"]
    G  = data.fields["Gbz"]
    κ  = abs.(dropdims(mean(P; dims=1); dims=1) ./ (2 .* (dropdims(mean(G; dims=1); dims=1) .+ ϵκ)))
    κ  = [mean(@view κ[max(1,k-4):min(end,k+4), max(1,n-6):min(end,n+6)]) for k in axes(κ, 1), n in axes(κ, 2)]
    κᵗ = vec(sum((κ[:, 1:end-1] + κ[:, 2:end]) .* reshape(diff(t), 1, :) / 2; dims=2)) / (t[end] - t[1])
    κᶻ = vec(sum(κ .* data.Δz; dims=1)) / sum(data.Δz)
    push!(profiles, (; label=variant.label, t, z, κ, κᵗ, κᶻ))
    data = nothing
    GC.gc()
end

# AB2-SE is saved every 18 steps rather than every hour, so each difference uses the saved time closest to each
# output time of WRK3-SE.
reference = only(filter(profile -> profile.label == "WRK3-SE", profiles))
nearest_records(t) = [argmin(abs.(t .- τ)) for τ in reference.t]
differences = [(; profile..., t=reference.t, Δκ=profile.κ[:, nearest_records(profile.t)] - reference.κ)
               for profile in profiles if profile.label in ("AB2-SE", "MRK4-SE", "WRK3-UP")]

plot_internal_tide(profiles, differences) #src
# ![Internal-tide numerical diffusivity](internal_tide_diffusivity.png)
