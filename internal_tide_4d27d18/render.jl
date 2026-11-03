# CELL 0
using Pkg
push!(LOAD_PATH, joinpath(dirname(Base.active_project()), "examples"))
using Suppressor
@suppress begin
    using Oceananigans, CairoMakie, TimestepperTestCases, Statistics
    using LaTeXStrings, Statistics
end

# CELL 1
# ============================================================================================
# The six time-stepping variants from internal_tide/test.jl, in canonical order.
# `free_surface` is the runner's `free_surface_name` (AB2 uses the default "split_free_surface").
# File loaded: folder * "internal_tide_" * timestepper * "_" * free_surface * ".jld2".
# Update `folder` to wherever the runs wrote their output (the runs are pending).
# ============================================================================================
folder = @__DIR__ * "/"
using JLD2
files = filter(name -> endswith(name, ".jld2"), readdir(folder; join=true))
@assert length(files) == 13 "Expected thirteen completed internal-tide simulations"
for file in files
    jldopen(file) do output
        times = output["timeseries/t"]
        @assert maximum(times[key] for key in keys(times)) ≈ 40 * 86400 "Incomplete simulation: $file"
    end
end

# Derived from the package so the notebook cannot drift from the runners.
variants = [(; label = d.label, timestepper = string(d.timestepper), free_surface = d.label) for d in discretizations()]

labels     = [v.label for v in variants]
styles     = [plot_style(v.label) for v in variants]
colors     = [s.color     for s in styles]
linestyles = [s.linestyle for s in styles]
cases  = [@suppress TimestepperTestCases.load_internal_tide(folder, v.timestepper, v.free_surface) for v in variants]

# CELL 2
using LaTeXStrings

ticks(vec) = (vec, latexstring.(string.(vec)))

function running_mean(v, points)
    n  = length(v)
    rm = zeros(length(v) - 2points+1)
    for i in points+1:n-points
        rm[i-points] = mean(v[i - points:i+points])
    end
    return rm[1:end-1]
end

path = mkpath(joinpath(@__DIR__, "figures")) * "/"

# CELL 3
# end-state fields for the snapshot figure, one entry per variant
uprime = [Field(c[:u][end] / c[:VFCC][end] - mean(c[:u][end] / c[:VFCC][end])) for c in cases]
ηend   = [c[:η][end] for c in cases]
N2end  = [Field(∂z(c[:b][end] / c[:VCCC][end])) for c in cases]

bu = Float64.(interior(cases[1][:u][end]) .== 0); bu[bu .== 1] .= NaN
bw = Float64.(interior(N2end[1]) .== 0);          bw[bw .== 1] .= NaN
nothing


# CELL 4
xu, yu, zu = nodes(cases[1][:u])
xw, yw, zw = nodes(cases[1][:w])
xη, yη, zη = nodes(cases[1][:η])

# Only the cases that differ. WRK3-SE-SM05 is indistinguishable from WRK3-SE, which Figure diffusivity-hovmoller
# establishes quantitatively.
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
    Colorbar(ga[r+5:r+7, ncols+1], hmN, label = L"\text{N}^{2}\text{ } 10^{-4}\text{[s}^{-2}\text{]}", ticks = ([0.95, 0.00, 1.05] .* 1e-4, latexstring.(string.([0.95, 0.00, 1.05]))))
end

colgap!(ga, 10); rowgap!(ga, 10)
for row in 1:nrows-1; rowgap!(ga, 7row, 30); end
nothing
CairoMakie.save(path * "figure-internal-tide.pdf", fig)

# CELL 5
mvv = 8
zC  = cases[1][:Abx].grid.z.cᵃᵃᶜ[2:128] ./ 1000
fig = Figure(size = (950, 350))

xticks_list = [
    ([0, 1e-5, 2e-5],         [L"0", L"10^{-5}", L"2 \cdot 10^{-5}"]),
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
              title  = L"t = \text{%$(Oceananigans.Utils.prettytime(cases[1][:Abz].times[it]))}",
              ylabel = col == 1 ? L"\text{Depth [km]}" : "",
              yticks = col == 1 ? ticks([-2, -1, 0]) : ([-2, -1, 0], ["", "", ""]))
    for (i, c) in enumerate(cases[1:end])
        κ = κprofile(c, it)
        lines!(ax, κ, zC, color = (colors[i], 0.1))
        lines!(ax, running_mean(κ, mvv), running_mean(zC, mvv), color = colors[i], linewidth = 1.5, label = labels[i], linestyle = linestyles[i])
    end
    ylims!(ax, -2, 0)
    xlims!(ax, 0, 5e-5)
    if col == 1; global axref = ax; end
end
Legend(fig[1, 5], axref)
display(current_figure())
CairoMakie.save(path * "diffusivity-internal-tide.pdf", fig)


# CELL 7
difference = true
ϵκ         = 1e-3 * 1e-4
nz_win     = 4
nt_win     = 6
ncols      = 3
κfloor     = 1e-9
κmax       = 1.5e-5      # colour saturation; raise to desaturate further

itr   = 1:length(cases[1][:Abz].times)
tdays = cases[1][:Abz].times[itr] ./ 86400

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

κs = [boxmean(reduce(hcat, [κprof(c, it) for it in itr]), nz_win, nt_win) for c in cases]

iref   = findfirst(v -> v.label == "WRK3-SE", variants)
panels = difference ? [i for i in eachindex(cases) if i != iref] : collect(eachindex(cases))
fields = difference ? [κs[i] .- κs[iref] for i in panels] : [κs[i] for i in panels]
crange = difference ? (-κmax, κmax) : (0, 4e-5)
cmap   = difference ? Reverse(:RdBu) : :magma

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
           title = L"\text{Time mean: 0--40 days}", titlealign = :left, titlesize = 13)

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
              ylabel = col == 1 ? L"\text{Depth [km]}" : "",
              yticks = col == 1 ? zticks : ([-2, -1, 0], ["", "", ""]),
              xticks = tticks,
              xticksvisible = bottom, xticklabelsvisible = bottom,
              xlabel = bottom ? L"\text{time [days]}" : "")

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
display(current_figure())
CairoMakie.save(path * "diffusivity-hovmoller-internal-tide.pdf", fig)

# CELL 8
fig = Figure(size = (820, 300))
ax1 = Axis(fig[1, 1], title=L"\text{TKE}",                    xlabel=L"\text{time [days]}", xticks=ticks([0, 10, 20, 30, 40]), yticks=ticks([0.00, 0.02, 0.04]))
ax2 = Axis(fig[1, 2], title=L"\text{RPE}-\text{RPE}|_{t=0}",  xlabel=L"\text{time [days]}", xticks=ticks([0, 10, 20, 30, 40]), yticks=ticks([0.00, 0.2, 0.4, 0.6]))
ax3 = Axis(fig[1, 3], title=L"\text{APE}-\text{APE}|_{t=0}",  xlabel=L"\text{time [days]}", xticks=ticks([0, 10, 20, 30, 40]), yticks=ticks([0.00, -0.2, -0.4, -0.6, -0.8]))
ax4 = Axis(fig[1, 4], title=L"\eta^2 \text{ [cm}^2\text{]}",  xlabel=L"\text{time [days]}", xticks=ticks([0, 10, 20, 30, 40]), yticks=ticks([0, 10, 20, 30, 40]))
times = cases[1][:Abx].times ./ 86400

for (i, c) in enumerate(cases)
    # lines!(ax1, times, c[:KE] .- c[:MKE], color = (colors[i], 0.05))
    lines!(ax1, running_mean(times, 12), running_mean(c[:KE] .- c[:MKE], 12), color = colors[i], label = labels[i], linestyle = linestyles[i])
    lines!(ax2, running_mean(times[3:end], 12), running_mean(c[:RPE][3:end] .- c[:RPE][3], 12), color = colors[i])
    # lines!(ax3, times[3:end], -(c[:APE][3:end] .- c[:APE][3]), color = (colors[i], 0.05))
    lines!(ax3, running_mean(times[3:end], 12), -running_mean(c[:APE][3:end] .- c[:APE][3], 12), color = colors[i], label = labels[i], linestyle = linestyles[i])
    # lines!(ax4, times, c[:η2] .* 1e4, color = (colors[i], 0.05))
    lines!(ax4, running_mean(times, 12), running_mean(c[:η2] .* 1e4, 12), color = colors[i], label = labels[i], linestyle = linestyles[i])
end

xlims!(ax1, 0, 40); xlims!(ax2, 0, 40); xlims!(ax3, 0, 40); xlims!(ax4, 0, 40)
Legend(fig[1, 5], ax3)
nothing
CairoMakie.save(path * "energetics-internal-tide.pdf", fig)


# CELL 11
# ══ Temporal-resolution test: the rectified bottom flow (section 5.1) ═════════════════════════════
# WRK3-SE and SRK3-SE at the full step, Δt/4 and Δt/8, with the substep count reduced in proportion so the
# barotropic Courant number is unchanged. Colour is the composition, dash is the reduced step: solid and
# dashed lying on top of each other means converged in time. The averaging kernel keeps its moments
# exactly at 64, 16 and 8 substeps, so a degraded filter is not confounding the comparison.

refined = [(; label = "WRK3-SE", timestepper = "SplitRungeKutta3",    tag = "dt4"),
           (; label = "SRK3-SE", timestepper = "SSPRungeKutta3",      tag = "dt4"),
           (; label = "MRK4-SE", timestepper = "ModifiedRungeKutta4", tag = "dt4"),
           (; label = "WRK3-SE", timestepper = "SplitRungeKutta3",    tag = "dt8"),
           (; label = "SRK3-SE", timestepper = "SSPRungeKutta3",      tag = "dt8"),
           (; label = "MRK4-SE", timestepper = "ModifiedRungeKutta4", tag = "dt8")]

cases_refined = [@suppress TimestepperTestCases.load_internal_tide(folder, v.timestepper, v.label * "-" * v.tag) for v in refined]

"Time-mean over the last `days` of the x-averaged u', in m/s: the tide cancels, the rectified flow survives."
function rectified_profile(c; days = 10)
    t  = c[:u].times ./ 86400
    i0 = findfirst(≥(t[end] - days), t)
    acc = zeros(size(interior(c[:u][1]), 3)); n = 0
    for it in i0:length(t)
        u = interior(c[:u][it]) ./ interior(c[:VFCC][it])     # z-star: volumes vary in time
        acc .+= vec(mean(u .- mean(u), dims = 1)); n += 1
    end
    return acc ./ n
end

# CELL 12
# ══ Temporal-resolution test: rectified flow and irreversible mixing ══════════════════════════════
zu_km = zu ./ 1000
iW = findfirst(v -> v.label == "WRK3-SE", variants)
iS = findfirst(v -> v.label == "SRK3-SE", variants)
iM = findfirst(v -> v.label == "MRK4-SE", variants)

curves = [("WRK3-SE",       cases[iW],        :solid, :blue),
          ("WRK3-SE, Δt/4", cases_refined[1], :dash,  :blue),
          ("WRK3-SE, Δt/8", cases_refined[4], :dot,   :blue),
          ("SRK3-SE",       cases[iS],        :solid, :red),
          ("SRK3-SE, Δt/4", cases_refined[2], :dash,  :red),
          ("SRK3-SE, Δt/8", cases_refined[5], :dot,   :red),
          ("MRK4-SE",       cases[iM],        :solid, :green),
          ("MRK4-SE, Δt/4", cases_refined[3], :dash,  :green),
          ("MRK4-SE, Δt/8", cases_refined[6], :dot,   :green)]

"Time-mean of u' over x and t: the tide cancels, the rectified flow survives."
function rectified_profile(c; days = 10)
    t  = c[:u].times ./ 86400
    i0 = findfirst(≥(t[end] - days), t)
    acc = zeros(size(interior(c[:u][1]), 3)); n = 0
    for it in i0:length(t)
        u = interior(c[:u][it]) ./ interior(c[:VFCC][it])
        acc .+= vec(mean(u .- mean(u), dims = 1)); n += 1
    end
    return acc ./ n
end

uticks = ([0.0, 1e-3, 2e-3, 3e-3], [L"0", L"1\cdot 10^{-3}", L"2\cdot 10^{-3}", L"3\cdot 10^{-3}"])
tticks = ticks([0, 10, 20, 30, 40])

fig = Figure(size = (900, 300), fontsize = 14)

axu = Axis(fig[1, 1]; xlabel = L"\langle u^\prime \rangle_{x,t} \text{ [m s}^{-1}\text{]}",
           ylabel = L"\text{Depth [km]}", xticks = uticks,
           yticks = ticks([-2.0, -1.8, -1.6, -1.4, -1.2]))

for (lab, c, ls, clr) in curves
    lines!(axu, rectified_profile(c), zu_km; color = clr, linestyle = ls, linewidth = 2, label = latexstring(replace(lab, "Δt" => raw"$\Delta t$")))
end
vlines!(axu, 0; color = (:black, 0.3), linewidth = 0.8)
ylims!(axu, -2, -1.2)

axr = Axis(fig[1, 2]; xlabel = L"\text{time [days]}", xticks = tticks,
           ylabel = L"\text{RPE} - \text{RPE}|_{t=0}", yticks = ticks([-0.1, 0.0, 0.1, 0.2, 0.3]))

for (lab, c, ls, clr) in curves
    t = c[:RPE] isa AbstractVector ? (c[:Abx].times ./ 86400) : nothing
    lines!(axr, running_mean(t[3:end], 12), running_mean(c[:RPE][3:end] .- c[:RPE][3], 12);
           color = clr, linestyle = ls, linewidth = 2, label = latexstring(replace(lab, "Δt" => raw"$\Delta t$")))
end
xlims!(axr, 0, 40)
Legend(fig[1, 3], axr; framevisible = true, labelsize = 11)

display(current_figure())
CairoMakie.save(path * "rectified-internal-tide.pdf", fig)

# CELL 13
