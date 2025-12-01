using TimestepperTestCases, Oceananigans, GLMakie
using LaTeXStrings, Statistics

rl = TimestepperTestCases.load_internal_tide("internal_tide/", "QuasiAdamsBashforth2", "split_free_surface")
ri = TimestepperTestCases.load_internal_tide("internal_tide/", "SplitRungeKutta3", "implicit_free_surface")
al = TimestepperTestCases.load_internal_tide("internal_tide/", "SplitRungeKutta3", "split_free_surface")

url = Field(rl[:u][end] / rl[:VFCC][end] - mean(rl[:u][end] / rl[:VFCC][end]))
uri = Field(ri[:u][end] / ri[:VFCC][end] - mean(ri[:u][end] / ri[:VFCC][end]))
ual = Field(al[:u][end] / al[:VFCC][end] - mean(al[:u][end] / al[:VFCC][end]))

wrl = Field(rl[:w][end] / rl[:VCCF][end])
wri = Field(ri[:w][end] / rl[:VCCF][end])
wal = Field(al[:w][end] / rl[:VCCF][end])

bu = Float64.(interior(rl[:u][end]) .== 0)
bu[bu .== 1] .= NaN

bw = Float64.(interior(wrl) .== 0)
bw[bw .== 1] .= NaN

ηrl = rl[:η][end]
ηri = ri[:η][end]
ηal = al[:η][end]

Nrl = Field(∂z(rl[:b][end] / rl[:VCCC][end]))
Nri = Field(∂z(ri[:b][end] / ri[:VCCC][end]))
Nal = Field(∂z(al[:b][end] / al[:VCCC][end]))

xu, yu, zu = nodes(rl[:u])
xw, yw, zw = nodes(rl[:w])
xη, yη, zη = nodes(rl[:η])

xticks  = ([-800, -400, 0, 400, 800] .* 1000, latexstring.(string.([-800, -400, 0, 400, 800])))
xnticks = ([-800, -400, 0, 400, 800] .* 1000, ["", "", "", "", ""])

yticks = ([-2, -1.5, -1, -0.5, -0] .* 1000, latexstring.(string.([-2.0, -1.5, -1.0, -0.5, -0.0])))
ynticks = ([-2, -1.5, -1, -0.5, -0] .* 1000, ["", "", "", "", ""])

yeticks =  ([-1, -0, 1] .* 1e-1, latexstring.(string.([-10, 0, 10])))
yenticks = ([-1, -0, 1] .* 1e-1, ["", "", ""])


fig  = Figure(resolution = (1100, 450), fontsize = 18)
ga   = GridLayout(fig[1, 1])
ax01 = Axis(ga[1, 1],    title = L"\text{\textbf{QAB2-IM}}", xlabel ="", ylabel = L"\eta \text{ [cm]}", xticks=xnticks, yticks=yeticks)
ax02 = Axis(ga[1, 2],    title = L"\text{\textbf{RK3-IM}}", xlabel ="", ylabel = "", xticks=xnticks,   yticks=yenticks)
ax03 = Axis(ga[1, 3],    title = L"\text{\textbf{RK3-SE}}", xlabel ="", ylabel = "", xticks=xnticks,   yticks=yenticks)
ax11 = Axis(ga[2:4, 1],  xlabel ="", ylabel = L"\text{Depth [km]}", xticks=xnticks, yticks=yticks)
ax12 = Axis(ga[2:4, 2],  xlabel ="", ylabel = "",                   xticks=xnticks, yticks=ynticks)
ax13 = Axis(ga[2:4, 3],  xlabel ="", ylabel = "",                   xticks=xnticks, yticks=ynticks)
ax21 = Axis(ga[5:7, 1], xlabel =L"\text{x [km]}", ylabel = L"\text{Depth [km]}", xticks=xticks, yticks=yticks)
ax22 = Axis(ga[5:7, 2], xlabel =L"\text{x [km]}", ylabel ="",                    xticks=xticks, yticks=ynticks)
ax23 = Axis(ga[5:7, 3], xlabel =L"\text{x [km]}", ylabel ="",                    xticks=xticks, yticks=ynticks)

lines!(ax01, xη, interior(ηal, :, 1, 1), color = :grey, linewidth = 2)
lines!(ax02, xη, interior(ηri, :, 1, 1), color = :grey, linewidth = 2)
lines!(ax03, xη, interior(ηrl, :, 1, 1), color = :grey, linewidth = 2)

xlims!(ax01, -1e6, 1e6)
xlims!(ax02, -1e6, 1e6)
xlims!(ax03, -1e6, 1e6)
ylims!(ax01, -0.15, 0.15)
ylims!(ax02, -0.15, 0.15)
ylims!(ax03, -0.15, 0.15)

hm = heatmap!(ax11, xu, zu, interior(ual, :, 1, :) .+ bu[:, 1, :]; colormap=:bwr, colorrange = (-0.35, 0.35), nan_color=:black)
hm = heatmap!(ax12, xu, zu, interior(uri, :, 1, :) .+ bu[:, 1, :]; colormap=:bwr, colorrange = (-0.35, 0.35), nan_color=:black)
hm = heatmap!(ax13, xu, zu, interior(url, :, 1, :) .+ bu[:, 1, :]; colormap=:bwr, colorrange = (-0.35, 0.35), nan_color=:black)
Colorbar(ga[2:4, 4], hm, label = L"\text{u' [m s}^{-1}\text{]}", ticks = ([-0.2, 0.0, 0.2], latexstring.(string.([-0.2, 0.0, 0.2]))))
hm = heatmap!(ax21, xw, zw[2:127], interior(Nal, :, 1, 2:127) .+ bw[:, 1, 2:127]; colormap=:ice, colorrange = (0.00009, 0.00011), nan_color = :green)
hm = heatmap!(ax22, xw, zw[2:127], interior(Nri, :, 1, 2:127) .+ bw[:, 1, 2:127]; colormap=:ice, colorrange = (0.00009, 0.00011), nan_color = :green)
hm = heatmap!(ax23, xw, zw[2:127], interior(Nrl, :, 1, 2:127) .+ bw[:, 1, 2:127]; colormap=:ice, colorrange = (0.00009, 0.00011), nan_color = :green)
Colorbar(ga[5:7, 4], hm, label = L"\text{N}^{2}\text{ } 10^{-4}\text{[s}^{-2}\text{]}", ticks = ([0.95, 0.00, 1.05] .* 1e-4, latexstring.(string.([0.95, 0.00, 1.05]))))

colgap!(ga, 10)
rowgap!(ga, 10)