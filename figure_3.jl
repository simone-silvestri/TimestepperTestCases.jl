using TimestepperTestCases, Oceananigans, GLMakie

ASxr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "ASx")
ASyr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "ASy")
ASzr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "ASz")
ASxa = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "ASx")
ASya = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "ASy")
ASza = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "ASz")

ATxr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "ATx")
ATyr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "ATy")
ATzr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "ATz")
ATxa = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "ATx")
ATya = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "ATy")
ATza = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "ATz")

aSxa = [sum(ASxa[i]) for i in 1:length(ASxa)]
aSya = [sum(ASya[i]) for i in 1:length(ASya)]
aSza = [sum(ASza[i]) for i in 1:length(ASza)]
aSxr = [sum(ASxr[i]) for i in 1:length(ASxr)]
aSyr = [sum(ASyr[i]) for i in 1:length(ASyr)]
aSzr = [sum(ASzr[i]) for i in 1:length(ASzr)]

tSa = aSxa .+ aSya .+ aSza
tSr = aSxr .+ aSyr .+ aSzr

GSxr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "GSx")
GSyr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "GSy")
GSzr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "GSz")
GSxa = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "GSx")
GSya = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "GSy")
GSza = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "GSz")

GTxr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "GTx")
GTyr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "GTy")
GTzr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "GTz")
GTxa = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "GTx")
GTya = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "GTy")
GTza = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "GTz")

gSxa = [sum(GSxa[i]) for i in 1:length(GSxa)]
gSya = [sum(GSya[i]) for i in 1:length(GSya)]
gSza = [sum(GSza[i]) for i in 1:length(GSza)]
gSzr = [sum(GSzr[i]) for i in 1:length(GSzr)]
gSyr = [sum(GSyr[i]) for i in 1:length(GSyr)]
gSxr = [sum(GSxr[i]) for i in 1:length(GSxr)]

DSza = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "DSz")
DSzr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "DSz")

dSza = [sum(DSza[i]) for i in 1:length(DSza)]
dSzr = [sum(DSzr[i]) for i in 1:length(DSzr)]

ur = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "u")
ua = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "u")
vr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "v")
va = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "v")
wr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "w")
wa = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "w")

KEa = [sum(ua[i]^2) + sum(va[i]^2) + sum(wa[i]^2) for i in 1:length(ua)]
KEr = [sum(ur[i]^2) + sum(vr[i]^2) + sum(wr[i]^2) for i in 1:length(ur)]

using Oceananigans.BuoyancyFormulations: buoyancy

Ta = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "T")
Tr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "T")
Sa = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_QuasiAdamsBashforth2_CATKE.jld2", "S")
Sr = FieldTimeSeries("idealized_coast/idealized_coast_split_free_surface_averages_SplitRungeKutta3_CATKE.jld2", "S")

V = sum(Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Center), Oceananigans.Operators.volume, Ta.grid))

b  = SeawaterBuoyancy(equation_of_state=LinearEquationOfState())
ba = FieldTimeSeries{Nothing, Center, Center}(Ta.grid, Ta.times)
br = FieldTimeSeries{Nothing, Center, Center}(Tr.grid, Tr.times)

for i in 1:length(Ta)
    set!(ba[i], buoyancy(b, Ta.grid, (T = Ta[i] / V1, S = Sa[i] / V1)))
    set!(br[i], buoyancy(b, Tr.grid, (T = Tr[i] / V1, S = Sr[i] / V1)))
end

fig = Figure(size = (1000, 600))

u_v = @lift(interior(ur[$iter], 1, :, :))
v_v = @lift(interior(vr[$iter], 1, :, :))
w_v = @lift(interior(wr[$iter], 1, :, :))
b_v = @lift(interior(br[$iter], 1, :, :))

ax = Axis(fig[1, 1]) 
heatmap!(ax, y, z,  u_v, colorrange=(-80000, 80000))
contour!(ax, y, z,  b_v, levels=range(-0.24, -0.1748, length=20), color=:black, linewidth=1)
ax = Axis(fig[1, 2]) 
heatmap!(ax, y, z,  v_v, colorrange=(-60000, 60000))
contour!(ax, y, z,  b_v, levels=range(-0.24, -0.1748, length=20), color=:black, linewidth=1)
ax = Axis(fig[1, 3]) 
heatmap!(ax, y, zF, w_v, colorrange=(-200, 200))
contour!(ax, y, z,  b_v, levels=range(-0.24, -0.1748, length=20), color=:black, linewidth=1)

ASa = FieldTimeSeries{Nothing, Center, Center}(Ta.grid, Ta.times)
ASr = FieldTimeSeries{Nothing, Center, Center}(Tr.grid, Tr.times)
GSa = FieldTimeSeries{Nothing, Center, Center}(Ta.grid, Ta.times)
GSr = FieldTimeSeries{Nothing, Center, Center}(Tr.grid, Tr.times)

κba = FieldTimeSeries{Nothing, Center, Center}(Ta.grid, Ta.times)
κbr = FieldTimeSeries{Nothing, Center, Center}(Tr.grid, Tr.times)

for t in 1:length(Ta)
    set!(Aba[t], (ASxa[t] + ASya[t] + ASza[t]))
    set!(Abr[t], (ASxr[t] + ASyr[t] + ASzr[t]))
    set!(Gba[t], (GSxa[t] + GSya[t] + GSza[t]))
    set!(Gbr[t], (GSxr[t] + GSyr[t] + GSzr[t]))
end

ka = mean(ASza[600], dims = (1, 2))[1,1,:] ./ mean(GSza[600], dims=(1, 2))[1,1,:] ./ 2 ./ (721 - 600)
kr = mean(ASzr[600], dims = (1, 2))[1,1,:] ./ mean(GSzr[600], dims=(1, 2))[1,1,:] ./ 2 ./ (721 - 600)

for t in 600:721
    ka .+= mean(ASza[t], dims = (1, 2))[1,1,:] ./ mean(GSza[t], dims=(1, 2))[1,1,:] ./ 2 ./ (721 - 600)
    kr .+= mean(ASzr[t], dims = (1, 2))[1,1,:] ./ mean(GSzr[t], dims=(1, 2))[1,1,:] ./ 2 ./ (721 - 600)
end