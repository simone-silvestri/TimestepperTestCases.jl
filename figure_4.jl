using Oceananigans, TimestepperTestCases, GLMakie
using Oceananigans.Units
using LaTeXStrings, Statistics

rl = TimestepperTestCases.load_internal_tide("internal_tide/", "QuasiAdamsBashforth2", "split_free_surface")
ri = TimestepperTestCases.load_internal_tide("internal_tide/", "SplitRungeKutta3", "implicit_free_surface")
al = TimestepperTestCases.load_internal_tide("internal_tide/", "SplitRungeKutta3", "split_free_surface")

myrange = 3:length(rl[:u])

rlR = (rl[:RPE][myrange] .- rl[:RPE][1]) ./ rl[:RPE][1]
riR = (ri[:RPE][myrange] .- ri[:RPE][1]) ./ ri[:RPE][1]
alR = (al[:RPE][myrange] .- al[:RPE][1]) ./ al[:RPE][1]

rlE = (rl[:KE][myrange] .- rl[:MKE][myrange])
riE = (ri[:KE][myrange] .- ri[:MKE][myrange])
alE = (al[:KE][myrange] .- al[:MKE][myrange])

rlA = (rl[:APE][myrange] .- rl[:APE][1]) ./ rl[:APE][1]
riA = (ri[:APE][myrange] .- ri[:APE][1]) ./ ri[:APE][1]
alA = (al[:APE][myrange] .- al[:APE][1]) ./ al[:APE][1]

rlH = rl[:η2][myrange]
riH = ri[:η2][myrange]
alH = al[:η2][myrange]

# rlκ = rl[:aSt][myrange] ./ rl[:gSt][myrange] ./ 2
# riκ = ri[:aSt][myrange] ./ ri[:gSt][myrange] ./ 2
# alκ = al[:aSt][myrange] ./ al[:gSt][myrange] ./ 2

times = rl[:u].times[myrange] ./ 1days

# rlκ = interior(mean(rl[:κb][300], dims=1), 1, 1, :) ./ (400 - 300)
# riκ = interior(mean(ri[:κb][300], dims=1), 1, 1, :) ./ (400 - 300)
# alκ = interior(mean(al[:κb][300], dims=1), 1, 1, :) ./ (400 - 300)

# for t in 301:400
#     rlκ .+= interior(mean(rl[:κb][t], dims=1), 1, 1, :) ./ (400 - 300)
#     riκ .+= interior(mean(ri[:κb][t], dims=1), 1, 1, :) ./ (400 - 300)
#     alκ .+= interior(mean(al[:κb][t], dims=1), 1, 1, :) ./ (400 - 300)
# end

function running_mean(v, points)
    n  = length(v)
    rm = zeros(length(v) - 2points+1)
    for i in points+1:n-points
        rm[i-points] = mean(v[i - points:i+points])
    end
    return rm[1:end-1]
end

c1 = Makie.wong_colors()[1]
c2 = Makie.wong_colors()[2]
c3 = Makie.wong_colors()[5]

fig = Figure(resolution = (1200, 350), fontsize = 20)
ax  = Axis(fig[1, 1], 
           xlabel = L"\text{time [days]}", 
           title  = L"10^5 \frac{RPE - RPE(t=0)}{RPE(t=0)}\text{ [-]}", 
           xticks = ([0, 10, 20, 30, 40], latexstring.(string.([0, 10, 20, 30, 40]))),
           yticks = ([3, 3.25, 3.5, 3.75, 4].*1e-5, latexstring.(string.([3.00, 3.25, 3.50, 3.75, 4.00])))) 

lines!(ax, times, alR, color = (c1, 0.15), linestyle = :solid, linewidth = 2, label = L"\text{\textbf{QAB2}, split free surface}")
lines!(ax, times, riR, color = (c2, 0.15), linestyle = :solid, linewidth = 2, label = L"\text{\textbf{RK3}, implicit free surface}")
lines!(ax, times, rlR, color = (c3, 0.15), linestyle = :solid, linewidth = 2, label = L"\text{\textbf{RK3}, split free surface}")

lines!(ax, running_mean(times, 12), running_mean(alR, 12), color = c1, linewidth = 2, linestyle = :solid)
lines!(ax, running_mean(times, 12), running_mean(riR, 12), color = c2, linewidth = 2, linestyle = :solid)
lines!(ax, running_mean(times, 12), running_mean(rlR, 12), color = c3, linewidth = 2, linestyle = :solid)

xlims!(ax, 0, 40)

ax  = Axis(fig[1, 2], 
           xlabel = L"\text{time [days]}", 
           title = L"TKE", 
           xticks = ([0, 10, 20, 30, 40], latexstring.(string.([0, 10, 20, 30, 40]))),
           yticks = ([0.00, 0.01, 0.02, 0.03, 0.04] .* 0.1, latexstring.(string.([0.00, 0.001, 0.002, 0.003, 0.004])))) 
lines!(ax, times, alE, color = (c1, 0.15), linestyle = :solid, label = L"\text{\textbf{QAB2}, split free surface}")
lines!(ax, times, riE, color = (c2, 0.15), linestyle = :solid, label = L"\text{\textbf{RK3}, implicit free surface}")
lines!(ax, times, rlE, color = (c3, 0.15), linestyle = :solid, label = L"\text{\textbf{RK3}, split free surface}")

lines!(ax, running_mean(times, 12), running_mean(alE, 12), linewidth = 2, color = c1, linestyle = :solid)
lines!(ax, running_mean(times, 12), running_mean(riE, 12), linewidth = 2, color = c2, linestyle = :solid)
lines!(ax, running_mean(times, 12), running_mean(rlE, 12), linewidth = 2, color = c3, linestyle = :solid)
xlims!(ax, 0, 40)

ax  = Axis(fig[1, 3], 
           xlabel = L"\text{time [days]}", 
           title  = L"10^7 \frac{APE - APE(t=0)}{APE(t=0)}\text{ [-]}", 
           xticks = ([0, 10, 20, 30, 40], latexstring.(string.([0, 10, 20, 30, 40]))),
           yticks = ([0, 5, 10, 15, 20].*1e-5, latexstring.(string.([0, 5, 10, 15, 20]))))  

lines!(ax, times, alA, color = (c1, 0.15), linestyle = :solid)
lines!(ax, times, riA, color = (c2, 0.15), linestyle = :solid)
lines!(ax, times, rlA, color = (c3, 0.15), linestyle = :solid)

lines!(ax, running_mean(times, 12), running_mean(alA, 12), color = c1, linewidth = 2, linestyle = :solid, label = L"\text{\textbf{QAB2}, split free surface}")
lines!(ax, running_mean(times, 12), running_mean(riA, 12), color = c2, linewidth = 2, linestyle = :solid, label = L"\text{\textbf{RK3}, implicit free surface}")
lines!(ax, running_mean(times, 12), running_mean(rlA, 12), color = c3, linewidth = 2, linestyle = :solid, label = L"\text{\textbf{RK3}, split free surface}")
axislegend(ax, position=:rb, framevisible=false)
xlims!(ax, 0, 40)

# ax  = Axis(fig[1, 4], 
#            xlabel = L"\text{time [days]}", 
#            title  = L"10^5 \frac{APE - APE(t=0)}{APE(t=0)}\text{ [-]}", 
#            xticks = ([0, 10, 20, 30, 40], latexstring.(string.([0, 10, 20, 30, 40]))),
#            yticks = ([-2.5, 0, 2.5, 5, 7.5, 10].*1e-5, latexstring.(string.([-2.5, 0, 2.5, 5, 7.5, 10]))))  

# lines!(ax, times, alH, color = (c1, 0.15), linestyle = :solid, label = L"\text{\textbf{RK3}, split free surface}")
# lines!(ax, times, riH, color = (c2, 0.15), linestyle = :solid, label = L"\text{\textbf{QAB2}, split free surface}")
# lines!(ax, times, rlH, color = (c3, 0.15), linestyle = :solid, label = L"\text{\textbf{RK3}, implicit free surface}")

# lines!(ax, running_mean(times, 12), running_mean(alH, 12), color = c1, linewidth = 2, linestyle = :solid)
# lines!(ax, running_mean(times, 12), running_mean(riH, 12), color = c2, linewidth = 2, linestyle = :solid)
# lines!(ax, running_mean(times, 12), running_mean(rlH, 12), color = c3, linewidth = 2, linestyle = :solid)
# xlims!(ax, 0, 40)
