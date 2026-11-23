# # Equilibrated channel: numerical diffusivity and dissipation
#
# In this example, we draw the channel figures of the manuscript (Figures 14 and 15) from the `averages_<label>.jld2` files of
# the record, which hold the buoyancy-variance tendencies and squared buoyancy gradients averaged over 1800-day windows.
# Only the final window is used. Run it with `julia --project=examples examples/figures/channel.jl /path/to/channel`;
# without an argument, a missing file is downloaded from the record into `data/channel_simulation/`.

using TimestepperTestCases
using Oceananigans
using Statistics
include("../shared/data_sources.jl") #src
include("../shared/compact_data.jl") #src
include("../shared/plotting.jl") #src
include("../shared/channel_plots.jl") #src

directory = isempty(ARGS) ? nothing : first(ARGS)
labels = ["AB2-SE", "WRK3-SE-SM05", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM", "WRK3-UP"]

# ## Zonal means of the final window
#
# The budget terms live on different cell faces, so `prepare_channel` interpolates them to cell centers before averaging
# them zonally. The reduction is cached in `data/channel_simulation/`, and the log records which averaging window was read.

datasets = map(labels) do label
    data = compact_dataset("channel_simulation", "averages_$label.jld2", prepare_channel; directory)
    @info "Channel averaging window" case=label final_year=data.time / (360 * 86400)
    return data
end

# ## Numerical diffusivity
#
# At each depth, the numerical diffusivity is minus the variance tendency summed over latitude, divided by twice the summed
# squared gradient. The horizontal diffusivity uses the `x` and `y` terms only, the total one all three directions. The
# vertical diffusivity of the mixing closure follows in the same way from its variance dissipation `Dbz`, and the third
# panel shows its ratio to the total numerical diffusivity. This ratio is the quantity that matters in an equilibrated
# ocean: where it exceeds one, the parameterized physics sets the mixing; where it falls below one, the discretization
# does.

profiles = map(labels, datasets) do label, data
    Pₕ = data.fields["Abx"] + data.fields["Aby"]
    Pz = data.fields["Abz"]
    Gₕ = data.fields["Gbx"] + data.fields["Gby"]
    Gz = data.fields["Gbz"]
    κₕ = vec(-sum(Pₕ; dims=1) ./ (2sum(Gₕ; dims=1)))
    κ  = vec(-sum(Pₕ + Pz; dims=1) ./ (2sum(Gₕ + Gz; dims=1)))
    κᵖ = vec(-sum(data.fields["Dbz"]; dims=1) ./ (2sum(Gz; dims=1)))
    return (; label, data.z, κₕ, κ, ratio=κᵖ ./ κ, data.time)
end

plot_channel_diffusivity(profiles) #src
# ![Channel numerical diffusivity](channel_diffusivity.png)
#
# The horizontal numerical diffusivity is 1–10 m² s⁻¹ in every case, and the total one between 10⁻⁶ and 10⁻⁴ m² s⁻¹.
# WRK3-UP mixes the most, and its ratio stays at or below one: the low-order advection controls the mixing over the whole
# column. AB2-SE mixes more than WRK3-SE. In the other cases the ratio exceeds one, with peaks above ten in the interior,
# and the Runge–Kutta schemes reach the largest values.

# ## Variance-tendency sections
#
# The sections show the horizontal and vertical variance tendencies of WRK3-SE and, for AB2-SE, MRK4-SE and WRK3-UP,
# their differences from WRK3-SE. The differences keep the sign of the tendency: negative values mark more variance
# destruction than in WRK3-SE.

sections = map(labels, datasets) do label, data
    return (; label, data.y, data.z, Pₕ=data.fields["Abx"] + data.fields["Aby"], Pz=data.fields["Abz"])
end

reference = only(filter(section -> section.label == "WRK3-SE", sections))
differences = [(; section..., ΔPₕ=section.Pₕ - reference.Pₕ, ΔPz=section.Pz - reference.Pz) for section in sections if section.label in ("AB2-SE", "MRK4-SE", "WRK3-UP")]

plot_channel_dissipation(reference, differences) #src
# ![Channel variance-tendency differences](channel_dissipation.png)
#
# The horizontal terms dominate the variance destruction, while the vertical ones show patches of variance production.
# WRK3-UP destroys by far the most variance; AB2-SE destroys more than WRK3-SE in the horizontal, which explains its
# larger horizontal diffusivity, and MRK4-SE slightly less.
