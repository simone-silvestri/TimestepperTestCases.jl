# # Internal tide: manuscript figures
#
# Render snapshots, numerical diffusivity, energetics, and timestep-refinement comparisons from the fifteen
# completed 40-day simulations. Run with `julia --project=examples examples/figures/internal_tide_figures.jl internal_tide`.

using Oceananigans, CairoMakie, TimestepperTestCases, Statistics, LaTeXStrings, JLD2
include("../shared/plotting.jl") #src
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

labels     = [v.label for v in variants]
styles     = [plot_style(v.label) for v in variants]
colors     = [s.color     for s in styles]
linestyles = [s.linestyle for s in styles]
cases  = [TimestepperTestCases.load_internal_tide(folder, v.timestepper, v.free_surface) for v in variants]


ticks(vec) = (vec, latexstring.(string.(vec)))

function running_mean(v, points)
    n  = length(v)
    rm = zeros(length(v) - 2points+1)
    for i in points+1:n-points
        rm[i-points] = mean(v[i - points:i+points])
    end
    return rm[1:end-1]
end

path = mkpath(joinpath(folder, "figures")) * "/"

# end-state fields for the snapshot figure, one entry per variant
uprime = [Field(c[:u][end] / c[:VFCC][end] - mean(c[:u][end] / c[:VFCC][end])) for c in cases]
ηend   = [c[:η][end] for c in cases]
N2end  = [Field(∂z(c[:b][end] / c[:VCCC][end])) for c in cases]

bu = Float64.(interior(cases[1][:u][end]) .== 0); bu[bu .== 1] .= NaN
bw = Float64.(interior(N2end[1]) .== 0);          bw[bw .== 1] .= NaN
nothing


xu, yu, zu = nodes(cases[1][:u])
xw, yw, zw = nodes(cases[1][:w])
xη, yη, zη = nodes(cases[1][:η])

shown = [findfirst(v -> v.label == l, variants) for l in ("AB2-SE", "WRK3-SE", "SRK3-SE", "MRK4-SE", "WRK3-IM", "WRK3-UP")]
xhalf = maximum(xu)
xt  = [-400, 0, 400]
xticks   = (xt .* 1000, latexstring.(string.(xt)))
xnticks  = (xt .* 1000, fill("", length(xt)))
yticks   = ([-2, -1.5, -1, -0.5, -0] .* 1000, latexstring.(string.([-2.0, -1.5, -1.0, -0.5, -0.0])))
ynticks  = ([-2, -1.5, -1, -0.5, -0] .* 1000, ["", "", "", "", ""])
yeticks  = ([-1, -0, 1] .* 1e-1, latexstring.(string.([-10, 0, 10])))
yenticks = ([-1, -0, 1] .* 1e-1, ["", "", ""])

# every case is a block of seven layout rows: η on the first, u' on the next three and N² on the last three
ncols = 3
nrows = cld(length(shown), ncols)
fig   = Figure(size = (450 * ncols + 120, 520 * nrows), fontsize = 16)
ga    = GridLayout(fig[1, 1])

hmu = nothing
hmN = nothing
for (n, k) in enumerate(shown)
    row, j = cld(n, ncols), mod1(n, ncols)
    fst    = j == 1
    bottom = n + ncols > length(shown)
    r      = 7 * (row - 1)

    axη = Axis(ga[r+1, j],     title = variants[k].label, ylabel = fst ? L"\eta \text{ [cm]}"  : "", xticks = xnticks, yticks = fst ? yeticks : yenticks)
    axu = Axis(ga[r+2:r+4, j], ylabel = fst ? L"\text{Depth [km]}" : "", xticks = xnticks, yticks = fst ? yticks : ynticks)
    axN = Axis(ga[r+5:r+7, j], xlabel = bottom ? L"\text{x [km]}" : "", ylabel = fst ? L"\text{Depth [km]}" : "",
               xticks = bottom ? xticks : xnticks, yticks = fst ? yticks : ynticks)

    lines!(axη, xη, interior(ηend[k], :, 1, 1), color = :grey, linewidth = 2)
    ylims!(axη, -0.15, 0.15)

    global hmu = heatmap!(axu, xu, zu, interior(uprime[k], :, 1, :) .+ bu[:, 1, :]; colormap = :bwr, colorrange = (-0.35, 0.35), nan_color = :black)
    global hmN = heatmap!(axN, xw, zw[2:127], interior(N2end[k], :, 1, 2:127) .+ bw[:, 1, 2:127]; colormap = :ice, colorrange = (0.00009, 0.00011), nan_color = :green)

    for ax in (axη, axu, axN); xlims!(ax, -xhalf, xhalf); end
end

for row in 1:nrows
    r = 7 * (row - 1)
    Colorbar(ga[r+2:r+4, ncols+1], hmu, label = L"\text{u' [m s}^{-1}\text{]}", ticks = ([-0.2, 0.0, 0.2], latexstring.(string.([-0.2, 0.0, 0.2]))))
    Colorbar(ga[r+5:r+7, ncols+1], hmN, label = L"\text{N}^{2}\text{ } 10^{-4}\text{[s}^{-2}\text{]}", ticks = ([0.95, 1.00, 1.05] .* 1e-4, latexstring.(string.([0.95, 1.00, 1.05]))))
end

colgap!(ga, 10); rowgap!(ga, 10)
for row in 1:nrows-1; rowgap!(ga, 7row, 30); end
nothing
for extension in ("pdf", "png")
    CairoMakie.save(path * "figure-internal-tide." * extension, fig)
end
# ![Final surface height, velocity, and stratification](figure-internal-tide.png)

# AB2-SE is saved every 18 steps rather than every hour, so the runs are compared at the saved time closest to each
# hourly output of WRK3-SE
reference_times = cases[findfirst(v -> v.label == "WRK3-SE", variants)][:Abz].times
nearest_index(c, t) = argmin(abs.(c[:Abz].times .- t))

mvv = 8
zC  = cases[1][:Abx].grid.z.cᵃᵃᶜ[2:128] ./ 1000
fig = Figure(size = (950, 350))

xticks_list = [
    ([0, 2e-5, 4e-5],         [L"0", L"2\cdot 10^{-5}", L"4\cdot 10^{-5}"]),
    ([0, 2e-5, 4e-5, 6e-5],   [L"0", L"2\cdot 10^{-5}", L"4 \cdot 10^{-5}", L"6 \cdot 10^{-5}"]),
    ([0, 3e-5, 6e-5, 9e-5],   [L"0", L"3\cdot 10^{-5}", L"6 \cdot 10^{-5}", L"9 \cdot 10^{-5}"]),
    ([0, 5e-5, 1e-4, 1.5e-4], [L"0", L"5\cdot 10^{-5}", L"10^{-4}", L"1.5 \cdot 10^{-4}"]),
]
its = [100, 300, 600, 961]
ϵ   = 0

κprofile(c, it) = abs.(interior(Field(mean(c[:Abz][it] + c[:Abx][it], dims=1) / (mean(c[:Gbz][it], dims=1) + ϵ) / 2), 1, 1, 2:128))

axref = nothing
for (col, it) in enumerate(its)
    ax = Axis(fig[1, col];
              xlabel = L"\kappa \text{ [m}^2\text{s}^{-1}\text{]}",
              xticks = xticks_list[col],
              title  = L"t = \text{%$(Oceananigans.Utils.prettytime(reference_times[it]))}",
              ylabel = col == 1 ? L"\text{Depth [km]}" : "",
              yticks = col == 1 ? ticks([-2, -1, 0]) : ([-2, -1, 0], ["", "", ""]))
    for (i, c) in enumerate(cases[1:end])
        κ = κprofile(c, nearest_index(c, reference_times[it]))
        lines!(ax, κ, zC, color = (colors[i], 0.1))
        lines!(ax, running_mean(κ, mvv), running_mean(zC, mvv), color = colors[i], linewidth = 1.5, label = replace(labels[i], "SM05" => "F1"), linestyle = linestyles[i])
    end
    ylims!(ax, -2, 0)
    xlims!(ax, 0, 5e-5)
    if col == 1; global axref = ax; end
end
Legend(fig[1, 5], axref)

for extension in ("pdf", "png")
    CairoMakie.save(path * "diffusivity-internal-tide." * extension, fig)
end
# ![Numerical diffusivity profiles](diffusivity-internal-tide.png)


difference = true
ϵκ         = 1e-3 * 1e-4
vertical_window     = 4
temporal_window     = 6
ncols      = 3
κfloor     = 1e-9
κmax       = 1.5e-5

tdays = reference_times ./ 86400

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

κs = [boxmean(reduce(hcat, [κprof(c, nearest_index(c, t)) for t in reference_times]), vertical_window, temporal_window) for c in cases]

iref   = findfirst(v -> v.label == "WRK3-SE", variants)
panels = difference ? [i for i in eachindex(cases) if i != iref] : collect(eachindex(cases))
fields = difference ? [κs[i] .- κs[iref] for i in panels] : [κs[i] for i in panels]
crange = difference ? (-κmax, κmax) : (0, 4e-5)

tticks    = ticks([0, 10, 20, 30, 40])
zticks    = ticks([-2, -1, 0])
cbticks   = ([-κmax, -κmax/2, 0.0, κmax/2, κmax],
             [L"-1.5\cdot 10^{-5}", L"-7.5\cdot 10^{-6}", L"0", L"7.5\cdot 10^{-6}", L"1.5\cdot 10^{-5}"])
κticks    = ([0.0, 2e-5, 4e-5, 6e-5], [L"0", L"2\cdot 10^{-5}", L"4\cdot 10^{-5}", L"6\cdot 10^{-5}"])
κlogticks = ([1e-7, 1e-6, 1e-5, 1e-4], [L"10^{-7}", L"10^{-6}", L"10^{-5}", L"10^{-4}"])

nrows = cld(length(panels), ncols)
fig   = Figure(size = (1250, 640), fontsize = 15)

gtop = GridLayout(fig[1, 1])       # independent column sizing, so the log axis cannot
gbot = GridLayout(fig[2, 1])       # push the Hovmöller columns out of alignment

# ── summary row ───────────────────────────────────────────────────────────────────────────────────
averaging_durations = diff(tdays)
averaging_weights = [averaging_durations[1] / 2; (averaging_durations[1:end-1] .+ averaging_durations[2:end]) ./ 2; averaging_durations[end] / 2]
averaging_weights ./= sum(averaging_weights)
time_averaged_profiles = [values * averaging_weights for values in κs]

axg = Axis(gtop[1, 1]; xlabel = L"\text{time [days]}", xticks = tticks, yticks = κticks,
           ylabel = L"\langle \kappa \rangle_z \text{ [m}^2\text{s}^{-1}\text{]}")
axp = Axis(gtop[1, 2]; xlabel = L"\kappa \text{ [m}^2\text{s}^{-1}\text{]}",
           xscale = log10, xticks = κlogticks,
           ylabel = L"\text{Depth [km]}", yticks = zticks,
           title  = L"\text{Time mean: 0–40 days}", titlealign = :left, titlesize = 13)

for (i, c) in enumerate(cases)
    lines!(axg, tdays, vec(mean(κs[i], dims = 1));
           color = colors[i], linestyle = linestyles[i], linewidth = 1.8, label = latexstring(raw"\text{" * replace(labels[i], "-SM05" => "-F1") * "}"))
    lines!(axp, max.(time_averaged_profiles[i], κfloor), zC;
           color = colors[i], linestyle = linestyles[i], linewidth = 1.8)
end
xlims!(axg, 0, 40)
ylims!(axp, -2, 0); xlims!(axp, 1e-7, 2e-4)

Legend(gtop[1, 3], axg; framevisible = true, labelsize = 12, patchsize = (22, 12))
colsize!(gtop, 1, Relative(0.56))
colsize!(gtop, 2, Relative(0.24))

# ── Hovmöller panels ──────────────────────────────────────────────────────────────────────────────
hm = nothing
for (r, i) in enumerate(panels)
    row, col = cld(r, ncols), mod1(r, ncols)
    bottom   = r + ncols > length(panels)

    ax = Axis(gbot[row, col]; title = latexstring(raw"\text{" * replace(labels[i], "-SM05" => "-F1") * "}"), titlealign = :left, titlesize = 13, titlegap = 3,
              ylabel        = col == 1 ? L"\text{Depth [km]}" : "",
              yticks        = col == 1 ? zticks : ([-2, -1, 0], ["", "", ""]),
              xticks        = tticks,
              xticksvisible = bottom, xticklabelsvisible = bottom,
              xlabel        = bottom ? L"\text{time [days]}" : "")

    global hm = heatmap!(ax, tdays, zC, fields[r]'; colorrange = crange, colormap =  Reverse(:RdBu),
                         lowclip = :darkblue, highclip = :darkred, rasterize = 2)
    ylims!(ax, -2, 0); xlims!(ax, 0, 40)
end

Colorbar(gbot[1:nrows, ncols+1], hm; ticks = cbticks,
         label = difference ? L"\kappa - \kappa_{\text{WRK3-SE}} \text{ [m}^2\text{s}^{-1}\text{]}" :
                             L"\kappa \text{ [m}^2\text{s}^{-1}\text{]}")

for c in 1:ncols; colsize!(gbot, c, Auto(1.0)); end   # equal-width panels

rowsize!(fig.layout, 1, Relative(0.33))
rowgap!(gbot, 8); colgap!(gbot, 10)
rowgap!(fig.layout, 14)

for extension in ("pdf", "png")
    CairoMakie.save(path * "diffusivity-hovmoller-internal-tide." * extension, fig)
end
# ![Numerical diffusivity history and differences](diffusivity-hovmoller-internal-tide.png)

# compute_rpe_density integrates ρ z with ρ = ρ₀(1 - b/g): g/ρ₀ turns it into the buoyancy form -∫ b z, in m² s⁻²
g, ρ₀ = 9.80655, 1000.0
potential_energy_scale = 1e3 * g / ρ₀

fig = Figure(size = (820, 300))
ax1 = Axis(fig[1, 1], title=L"\text{TKE}\;[10^{-3}\,\mathrm{m^2\,s^{-2}}]",                    xlabel=L"\text{time [days]}", xticks=ticks([0, 10, 20, 30, 40]), yticks=WilkinsonTicks(4))
ax2 = Axis(fig[1, 2], title=L"\Delta\text{RPE}\;[10^{-3}\,\mathrm{m^2\,s^{-2}}]",  xlabel=L"\text{time [days]}", xticks=ticks([0, 10, 20, 30, 40]), yticks=ticks([0, 2, 4, 6]))
ax3 = Axis(fig[1, 3], title=L"\Delta\text{APE}\;[10^{-3}\,\mathrm{m^2\,s^{-2}}]",  xlabel=L"\text{time [days]}", xticks=ticks([0, 10, 20, 30, 40]), yticks=ticks([0, 2, 4, 6]))
ax4 = Axis(fig[1, 4], title=L"\eta^2 \text{ [cm}^2\text{]}",  xlabel=L"\text{time [days]}", xticks=ticks([0, 10, 20, 30, 40]), yticks=ticks([0, 10, 20, 30, 40]))
for (i, c) in enumerate(cases)
    times = c[:Abx].times ./ 86400
    lines!(ax1, running_mean(times, 12), 1e3 * running_mean(c[:KE] .- c[:MKE], 12), color = colors[i], label = replace(labels[i], "SM05" => "F1"), linestyle = linestyles[i])
    lines!(ax2, running_mean(times[3:end], 12), potential_energy_scale * running_mean(c[:RPE][3:end] .- c[:RPE][3], 12), color = colors[i], linestyle = linestyles[i])
    lines!(ax3, running_mean(times[3:end], 12), potential_energy_scale * running_mean(c[:APE][3:end] .- c[:APE][3], 12), color = colors[i], label = replace(labels[i], "SM05" => "F1"), linestyle = linestyles[i])
    lines!(ax4, running_mean(times, 12), running_mean(c[:η2] .* 1e4, 12), color = colors[i], label = replace(labels[i], "SM05" => "F1"), linestyle = linestyles[i])
end

xlims!(ax1, 0, 40); xlims!(ax2, 0, 40); xlims!(ax3, 0, 40); xlims!(ax4, 0, 40)
Legend(fig[1, 5], ax3)
nothing
for extension in ("pdf", "png")
    CairoMakie.save(path * "energetics-internal-tide." * extension, fig)
end
# ![Energy histories](energetics-internal-tide.png)

production_only && exit()


# ══ Temporal-resolution test: rectified flow, numerical diffusivity and TKE (section 5.1) ══════════════════════
#
# The curves at the nominal step are the production runs.

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

depths_in_kilometers = zu ./ 1000

"Time-mean of u' over x and t: the tide cancels, the rectified flow survives."
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

"Time mean of the box-averaged diffusivity profile, as in the summary row of the Hovmöller figure."
function time_mean_diffusivity(c)
    t = c[:Abz].times ./ 86400
    κ = boxmean(reduce(hcat, [κprof(c, it) for it in eachindex(t)]), vertical_window, temporal_window)
    Δ = diff(t)
    w = [Δ[1] / 2; (Δ[1:end-1] .+ Δ[2:end]) ./ 2; Δ[end] / 2]
    return κ * (w ./ sum(w))
end

refinement_label(lab) = latexstring(raw"\text{" * replace(lab, "-" => "‐", "Δt" => raw"}\Delta t\text{") * "}")

refinement = [(; label = refinement_label(lab), color = clr, linestyle = ls,
                 u = rectified_profile(c), κ = max.(time_mean_diffusivity(c), κfloor),
                 t = c[:Abx].times ./ 86400, tke = 1e3 .* (c[:KE] .- c[:MKE]))
              for (lab, c, ls, clr) in curves]

uticks = ([0.0, 1e-3, 2e-3, 3e-3], [L"0", L"1", L"2", L"3"])
tticks = ticks([0, 20, 40])

fig = Figure(size = (1600, 540), fontsize = 26)

axu = Axis(fig[1, 1]; xlabel = L"\langle u^\prime \rangle_{x,t} \text{ [10}^{-3}\text{ m s}^{-1}\text{]}", xticks = uticks,
           ylabel = L"\text{Depth [km]}", yticks = ticks([-2.0, -1.6, -1.2]))
axκ = Axis(fig[1, 2]; xlabel = L"\kappa \text{ [m}^2\text{s}^{-1}\text{]}", xscale = log10, xticks = κlogticks,
           ylabel = L"\text{Depth [km]}", yticks = zticks)
axk = Axis(fig[1, 3]; xlabel = L"\text{time [days]}", xticks = tticks,
           ylabel = L"\text{TKE}\;[10^{-3}\,\mathrm{m^2\,s^{-2}}]", yticks = WilkinsonTicks(3))

for r in refinement
    style = (; color = r.color, linestyle = r.linestyle, linewidth = 2)
    lines!(axu, r.u, depths_in_kilometers; style...)
    lines!(axκ, r.κ, zC; style...)
    lines!(axk, running_mean(r.t, 12), running_mean(r.tke, 12); label = r.label, style...)
end
vlines!(axu, 0; color = (:black, 0.3), linewidth = 0.8)
ylims!(axu, -2, -1.2)
ylims!(axκ, -2, 0); xlims!(axκ, 1e-7, 2e-4)
xlims!(axk, 0, 40)

Legend(fig[1, 4], axk; framevisible = true, labelsize = 24, patchsize = (40, 14), rowgap = 8, patchlabelgap = 10)


for extension in ("pdf", "png")
    CairoMakie.save(path * "rectified-internal-tide." * extension, fig)
end
# ![Timestep refinement](rectified-internal-tide.png)

