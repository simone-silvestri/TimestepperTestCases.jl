# # Internal tide: manuscript figures
#
# In this example, we draw the internal-tide figures of the manuscript (Figures 10–13) from the fifteen 40-day
# simulations archived on Zenodo: the seven configurations of Table 2 at their production time step, and AB2-SE, WRK3-SE, SRK3-SE and MRK4-SE
# at one fourth and one eighth of it. Run it with
# `julia --project=examples examples/figures/internal_tide_figures.jl /path/to/internal_tide`; the figures are written
# to `/path/to/internal_tide/figures/`. The option `--production-only` reads only the seven production runs and skips
# the timestep-refinement figure.

using Oceananigans, CairoMakie, TimestepperTestCases, Statistics, LaTeXStrings, JLD2
include("../shared/plotting.jl") #src
include("../shared/internal_tide_plots.jl") #src

# ## Loading the simulations
#
# We first check that every run reached day 40, then load each discretization. The loader returns the saved fields as
# `FieldTimeSeries`, together with the energy diagnostics computed from them.

folder = abspath(isempty(ARGS) ? joinpath(@__DIR__, "..", "..", "internal_tide") : first(ARGS)) * "/"
production_only = "--production-only" in ARGS
files = filter(name -> endswith(name, ".jld2"), readdir(folder; join=true))
production_only && filter!(name -> !occursin(r"-dt[48]\.jld2$", name), files)
@assert length(files) == (production_only ? 7 : 15) "Unexpected number of internal-tide simulations"
for file in files
    jldopen(file) do output
        times = output["timeseries/t"]
        @assert maximum(times[key] for key in keys(times)) > 40 * 86400 - 3600 "Incomplete simulation: $file"
    end
end

variants = [(; label = d.label, timestepper = string(d.timestepper), free_surface = d.label) for d in discretizations()]
cases    = [TimestepperTestCases.load_internal_tide(folder, v.timestepper, v.free_surface) for v in variants]

labels     = [v.label for v in variants]
styles     = [plot_style(v.label) for v in variants]
colors     = [s.color     for s in styles]
linestyles = [s.linestyle for s in styles]
path       = mkpath(joinpath(folder, "figures")) * "/"

# ## The flow after 40 days
#
# The saved velocity and buoyancy carry the volume of their cell, so they are divided by `VFCC` and `VCCC` before use.
# Our first figure shows the free surface, the velocity fluctuation ``u' = u - V^{-1} \int_V u \, \mathrm{d}V`` and the
# stratification ``N^2 = \partial_z b`` at the final output of six configurations.

uprime = [Field(c[:u][end] / c[:VFCC][end] - mean(c[:u][end] / c[:VFCC][end])) for c in cases]
ηend   = [c[:η][end] for c in cases]
N2end  = [Field(∂z(c[:b][end] / c[:VCCC][end])) for c in cases]

plot_internal_tide_snapshot(variants, cases, uprime, ηend, N2end, path) #src
# ![Final surface height, velocity, and stratification](figure-internal-tide.png)
#
# The Runge–Kutta solutions are more energetic than AB2-SE: both ``u'`` and ``N^2`` keep sharper structures, a first sign
# that the Runge–Kutta schemes damp the internal tide less.

# ## Numerical diffusivity profiles
#
# How much does the discretization mix? A diffusive flux ``-\kappa \nabla b`` destroys buoyancy variance at the rate
# ``2 \kappa |\nabla b|^2``, so dividing the variance destruction that the advection scheme actually produces by twice
# the squared gradient tells us the diffusivity the discretization applies,
#
# ```math
# \kappa = \left| \frac{\langle P_x + P_z \rangle_x}{2 \langle G_z \rangle_x} \right| ,
# ```
#
# where ``P_x`` and ``P_z`` are the saved variance tendencies `Abx` and `Abz`, ``G_z`` is the squared vertical gradient
# `Gbz`, and ``\langle \cdot \rangle_x`` is a horizontal mean. In this uniformly stratified flow the vertical gradient
# dominates the squared gradient.
# AB2-SE is saved every 18 steps rather than every hour, so every run is read at the saved time closest to each hourly
# output of WRK3-SE.

reference_times = cases[findfirst(v -> v.label == "WRK3-SE", variants)][:Abz].times
nearest_index(c, t) = argmin(abs.(c[:Abz].times .- t))
zC = cases[1][:Abx].grid.z.cᵃᵃᶜ[2:128] ./ 1000

κprofile(c, it) = abs.(interior(Field(mean(c[:Abz][it] + c[:Abx][it], dims=1) / mean(c[:Gbz][it], dims=1) / 2), 1, 1, 2:128))

profile_records = [100, 300, 600, 961]
profiles = [[κprofile(c, nearest_index(c, reference_times[it])) for c in cases] for it in profile_records]

plot_internal_tide_diffusivity_profiles(profiles, reference_times[profile_records], zC, labels, colors, linestyles, path) #src
# ![Numerical diffusivity profiles](diffusivity-internal-tide.png)
#
# AB2-SE is about twice as diffusive as WRK3-SE, and in the interior it mixes as much as the low-order spatial scheme of
# WRK3-UP, even though it takes a step 0.29 times as long. Near the bottom, where the staircase bathymetry lowers the
# order of the advection operators, all schemes mix alike.

# ## Diffusivity history
#
# For the depth–time picture the denominator is regularized by ``10^{-7}`` so that unstratified levels do not divide by
# zero, and each profile history is smoothed by a box mean over ±4 levels and ±6 outputs. The time mean of the profiles
# uses trapezoidal weights over the 40 days.

ϵκ = 1e-3 * 1e-4
κprof(c, it) = abs.(interior(Field(mean(c[:Abz][it] + c[:Abx][it], dims=1) /
                                   (mean(c[:Gbz][it], dims=1) + ϵκ) / 2), 1, 1, 2:128))

function boxmean(A, nz, nt)
    nzA, ntA = size(A)
    B = similar(A)
    for j in 1:ntA, i in 1:nzA
        B[i, j] = mean(@view A[max(1, i-nz):min(nzA, i+nz), max(1, j-nt):min(ntA, j+nt)])
    end
    return B
end

vertical_window = 4
temporal_window = 6
κs = [boxmean(reduce(hcat, [κprof(c, nearest_index(c, t)) for t in reference_times]), vertical_window, temporal_window) for c in cases]

tdays = reference_times ./ 86400
averaging_durations = diff(tdays)
averaging_weights = [averaging_durations[1] / 2; (averaging_durations[1:end-1] .+ averaging_durations[2:end]) ./ 2; averaging_durations[end] / 2]
averaging_weights ./= sum(averaging_weights)
time_averaged_profiles = [values * averaging_weights for values in κs]

# The lower panels show each configuration minus WRK3-SE.

reference_index = findfirst(v -> v.label == "WRK3-SE", variants)
panels = [i for i in eachindex(cases) if i != reference_index]
fields = [κs[i] .- κs[reference_index] for i in panels]

plot_internal_tide_diffusivity_history(tdays, zC, κs, time_averaged_profiles, panels, fields, labels, colors, linestyles, path) #src
# ![Numerical diffusivity history and differences](diffusivity-hovmoller-internal-tide.png)
#
# WRK3-SE-F1, SRK3-SE and MRK4-SE are indistinguishable from WRK3-SE: neither the barotropic treatment nor the
# arrangement of the Runge–Kutta stages leaves a mark on the mixing. The largest differences are those of AB2-SE and
# WRK3-UP.

# ## Energetics
#
# Numerical mixing and numerical damping also show up in the energy budget. The turbulent kinetic energy (TKE) is the
# kinetic energy minus that of the horizontally averaged flow. The reference potential energy (RPE) is the potential
# energy of the adiabatically sorted state; only irreversible mixing raises it, so its growth measures the spurious
# mixing accumulated by each discretization. The reference and available potential energies come from
# `compute_rpe_density`, which sorts the density anomaly ``\rho' = -\rho_0 (b - \bar b) / g`` into the basin and
# integrates ``\rho' z_\star``. Each cell enters the sort as 32 vertical slices carrying a limited linear profile: the
# initial buoyancy is uniform along each grid level, and whole cells would make every level a block of identical parcels
# whose reference energy drops, by up to ``N^2 \Delta z^2 / 12``, as soon as the tide heaves the isopycnals. The factor
# ``g / \rho_0`` turns the energies into the buoyancy form ``-\int (b - \bar b) z_\star``, in m² s⁻², shown as changes
# from the third output. The curves are running means over 25 outputs.

g, ρ₀ = 9.80655, 1000.0
potential_energy_scale = 1e3 * g / ρ₀

plot_internal_tide_energetics(cases, labels, colors, linestyles, potential_energy_scale, path) #src
# ![Energy histories](energetics-internal-tide.png)
#
# The TKE of AB2-SE saturates after about 20 days, well below that of the Runge–Kutta runs, and the RPE of AB2-SE grows
# about 1.7 times as much as theirs. Only WRK3-UP mixes more, as its higher TKE feeds the low-order advection.

# ## Timestep refinement
#
# AB2-SE, WRK3-SE, SRK3-SE and MRK4-SE also ran at one fourth and one eighth of their production step; the production
# runs are the curves at the nominal step. The figure compares the rectified flow, the time-mean numerical diffusivity
# and the turbulent kinetic energy.

production_only && exit()

refinement_schemes = [("AB2-SE",  "QuasiAdamsBashforth2", :black),
                      ("WRK3-SE", "SplitRungeKutta3",     :blue),
                      ("SRK3-SE", "SSPRungeKutta3",       :red),
                      ("MRK4-SE", "ModifiedRungeKutta4",  :green)]

function refinement_case(label, timestepper, n)
    n == 1 && return cases[findfirst(v -> v.label == label, variants)]
    return TimestepperTestCases.load_internal_tide(folder, timestepper, "$label-dt$n")
end

curves = [(n == 1 ? label : "$label, Δt/$n", refinement_case(label, timestepper, n), linestyle, color)
          for (label, timestepper, color) in refinement_schemes
          for (n, linestyle) in ((1, :solid), (4, :dash), (8, :dot))]

# The rectified flow is the mean of ``u'`` over ``x`` and over the last ten days, in which the tide cancels. The
# diffusivity is the trapezoidal time mean of the box-averaged profiles used above.

function rectified_profile(c; days = 10)
    t   = c[:u].times ./ 86400
    i0  = findfirst(≥(t[end] - days), t)
    acc = zeros(size(interior(c[:u][1]), 3)); n = 0
    for it in i0:length(t)
        u = interior(c[:u][it]) ./ interior(c[:VFCC][it])
        acc .+= vec(mean(u .- mean(u), dims = 1)); n += 1
    end
    return acc ./ n
end

function time_mean_diffusivity(c)
    t = c[:Abz].times ./ 86400
    κ = boxmean(reduce(hcat, [κprof(c, it) for it in eachindex(t)]), vertical_window, temporal_window)
    Δ = diff(t)
    w = [Δ[1] / 2; (Δ[1:end-1] .+ Δ[2:end]) ./ 2; Δ[end] / 2]
    return κ * (w ./ sum(w))
end

κfloor = 1e-9
refinement = [(; label, color, linestyle, u = rectified_profile(c), κ = max.(time_mean_diffusivity(c), κfloor),
                 t = c[:Abx].times ./ 86400, tke = 1e3 .* (c[:KE] .- c[:MKE]))
              for (label, c, linestyle, color) in curves]

depths_in_kilometers = nodes(cases[1][:u])[3] ./ 1000
plot_internal_tide_refinement(refinement, depths_in_kilometers, zC, path) #src
# ![Timestep refinement](rectified-internal-tide.png)
#
# The Runge–Kutta curves do not move when the step is refined, while AB2-SE approaches them as its step shrinks: its
# lower energy and larger mixing come from the time discretization. Even at one eighth of its step, 0.41 minutes,
# AB2-SE does not reach the Runge–Kutta solutions at 11.3 minutes.
