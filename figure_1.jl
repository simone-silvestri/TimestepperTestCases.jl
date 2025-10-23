qab2(k, Δt, e) = abs((1 - (1.5 + e) * im * k * Δt + sqrt(1 - (1 - 2e) * im * k * Δt - (3/2+e)^2 * (k * Δt)^2)) / 2)
f(k, Δt) = abs(1 - im * k * Δt + (im * k * Δt)^2 / 2 - (im * k * Δt)^3 / 6)
g(k, Δt) = abs(1 - im * k * Δt + (im * k * Δt)^2 / 2 - (im * k * Δt)^3 / 6 + (im * k * Δt)^4 / 24)
# h(k, Δt) = abs(1 - im * k * Δt + (im * k * Δt)^2 / 2 - (im * k * Δt)^3 / 6 + (im * k * Δt)^4 / 24 - (im * k * Δt)^5 / 120)
# l(k, Δt) = abs(1 - im * k * Δt + (im * k * Δt)^2 / 2 - (im * k * Δt)^3 / 6 + (im * k * Δt)^4 / 24 - (im * k * Δt)^5 / 120 + (im * k * Δt)^6 / 720)

using GLMakie
using LaTeXStrings

fig = Figure(size = (600, 300))
ax = Axis(fig[1:2, 1:2], 
          ylabel = L"\text{Amplification factor \alpha}", 
          xlabel = L"\text{CFL} = k\cdot \Delta t", 
          xticks = ([0, 0.25, 0.5, 0.75, 1], latexstring.(string.([0, 0.25, 0.5, 0.75, 1]))), 
          yticks = ([0.99, 0.995, 1.0, 1.005, 1.01], latexstring.(string.([0.99, 0.995, 1.0, 1.005, 1.01]))))

C = 0:0.01:3

lines!(C, qab2.(1, C, 0.05), color = :grey,  linestyle = :dash,  label = L"QAB2, \ \epsilon = 0.05")
lines!(C, qab2.(1, C, 0.1),  color = :grey,  linestyle = :solid, label = L"QAB2, \ \epsilon = 0.1")
lines!(C, f.(1, C),          color = :black, linestyle = :dash,  label = L"RK, \ m = 3")
lines!(C, g.(1, C),          color = :black, linestyle = :solid, label = L"RK, \ m = 4")

ylims!(ax, 0.985, 1.015)
xlims!(ax, 0, 3)

# ax2 = Axis(fig[1, 3], 
#           xticks = ([0, 1, 2.0, 3.0], latexstring.(string.([0, 1, 2, 3]))), 
#           yticks = ([0.9, 0.95, 1.0, 1.05, 1.1], latexstring.(string.([0.9, 0.95, 1.0, 1.05, 1.1]))))
# lines!(ax2, C .* 3, f.(1, C .* 3), color = :black, linestyle = :dash)
# lines!(ax2, C .* 3, g.(1, C .* 3), color = :black, linestyle = :solid)

# ylims!(ax2, 0.90, 1.1)
# xlims!(ax2, 0, 3)

Legend(fig[2, 3], ax)