using TimestepperTestCases, Oceananigans, GLMakie
using Oceananigans.AbstractOperations: grid_metric_operation
using Oceananigans.Operators: volume
using FFTW

function add_surface_dynamics!(case::Dict, folder, closure, suffix, timestepper)

    path_surf = folder * "idealized_coast_" * suffix * "surface_" * timestepper * "_" * closure * ".jld2"

    case[:uS] = FieldTimeSeries(path_surf, "u")
    case[:vS] = FieldTimeSeries(path_surf, "v")

    case[:specU] = []
    case[:specV] = []

    Nt = length(case[:uS])

    for t in 1:Nt
        usurf = Field(case[:uS][t] / grid_metric_operation((Face, Center, Face), volume, case[:uS].grid))
        vsurf = Field(case[:vS][t] / grid_metric_operation((Center, Face, Face), volume, case[:vS].grid))

        spec_u, freq = y_average_spectra(usurf, 1:size(usurf, 1), 1:size(usurf, 2); k=1)
        spec_v, freq = y_average_spectra(vsurf, 1:size(vsurf, 1), 1:size(vsurf, 2); k=1)

        push!(case[:specU], spec_u)
        push!(case[:specV], spec_v)
    end

    return case
end

function power_spectrum_1d_x(var, x)
    Nx  = length(x)
    Nfx = Int64(Nx)

    spectra = zeros(ComplexF64, Int(Nfx/2))
    dx = x[2] - x[1]

    freqs = fftfreq(Nfx, 1.0 / dx) # 0,+ve freq,-ve freqs (lowest to highest)
    freqs = freqs[1:Int(Nfx/2)] .* 2.0 .* π

    fourier      = fft(var) / Nfx
    spectra[1]  += fourier[1] .* conj(fourier[1])

    for m in 2:Int(Nfx/2)
        spectra[m] += 2.0 * fourier[m] * conj(fourier[m]) # factor 2 for neg freq contribution
    end

    return spectra, freqs
end

function y_average_spectra(var::Field, irange, jrange; k = 1, spectra = power_spectrum_1d_x)

    xdomain = xnodes(var)[irange]
    ydomain = ynodes(var)[jrange]
    spec, freq = spectra(interior(var, irange, jrange[1], k), xdomain)

    for j in jrange[2:end]
        spec .+= spectra(interior(var, irange, j, k), xdomain)[1]
    end

    spec .*= 1 / length(jrange)

    return spec, freq
end

function load_case(folder, closure, suffix, timestepper)
    path_avg = folder * "idealized_coast_" * suffix * "averages_" * timestepper * "_" * closure * ".jld2"

    case = Dict()
    case[:ASx] = FieldTimeSeries(path_avg, "ASx")
    case[:ASy] = FieldTimeSeries(path_avg, "ASy")
    case[:ASz] = FieldTimeSeries(path_avg, "ASz")
    case[:ATx] = FieldTimeSeries(path_avg, "ATx")
    case[:ATy] = FieldTimeSeries(path_avg, "ATy")
    case[:ATz] = FieldTimeSeries(path_avg, "ATz")

    Nt = length(case[:ASx])

    grid = case[:ASx].grid
    VFCC = grid_metric_operation((Face,   Center, Center), volume, grid)
    VCFC = grid_metric_operation((Center, Face,   Center), volume, grid)
    VCCF = grid_metric_operation((Center, Center, Face),   volume, grid)
    VCCC = grid_metric_operation((Center, Center, Center), volume, grid)

    case[:aSx] = [sum(case[:ASx][i] / VFCC) for i in 1:Nt] ./ prod(size(grid)) 
    case[:aSy] = [sum(case[:ASy][i] / VCFC) for i in 1:Nt] ./ prod(size(grid))
    case[:aSz] = [sum(case[:ASz][i] / VCCF) for i in 1:Nt] ./ prod(size(grid))
    case[:aSt] = case[:aSx] .+ case[:aSy] .+ case[:aSz]

    case[:GSx] = FieldTimeSeries(path_avg, "GSx")
    case[:GSy] = FieldTimeSeries(path_avg, "GSy")
    case[:GSz] = FieldTimeSeries(path_avg, "GSz")
    case[:GTx] = FieldTimeSeries(path_avg, "GTx")
    case[:GTy] = FieldTimeSeries(path_avg, "GTy")
    case[:GTz] = FieldTimeSeries(path_avg, "GTz")
    case[:DSz] = FieldTimeSeries(path_avg, "DSz")

    case[:gSx] = [sum(case[:GSx][i] / VFCC) for i in 1:Nt] ./ prod(size(grid))
    case[:gSy] = [sum(case[:GSy][i] / VCFC) for i in 1:Nt] ./ prod(size(grid))
    case[:gSz] = [sum(case[:GSz][i] / VCCF) for i in 1:Nt] ./ prod(size(grid))
    case[:dSz] = [sum(case[:DSz][i] / VCCC) for i in 1:Nt] ./ prod(size(grid))
    case[:gSt] = case[:gSx] .+ case[:gSy] .+ case[:gSz]

    case[:uX] = FieldTimeSeries(path_avg, "u")
    case[:vX] = FieldTimeSeries(path_avg, "v")
    case[:wX] = FieldTimeSeries(path_avg, "w")
    case[:TX] = FieldTimeSeries(path_avg, "T")
    case[:SX] = FieldTimeSeries(path_avg, "S")
    case[:bX] = FieldTimeSeries(path_avg, "b")

    # case[:u] = FieldTimeSeries(path_avg, "u", backend=OnDisk())
    # case[:v] = FieldTimeSeries(path_avg, "v", backend=OnDisk())
    # case[:w] = FieldTimeSeries(path_avg, "w", backend=OnDisk())
    # case[:T] = FieldTimeSeries(path_avg, "T", backend=OnDisk())
    # case[:S] = FieldTimeSeries(path_avg, "S", backend=OnDisk())
    # case[:b] = FieldTimeSeries(path_avg, "b", backend=OnDisk())

    # u2(i) = (case[:u][i] / VFCC)^2 * VFCC
    # v2(i) = (case[:v][i] / VCFC)^2 * VCFC
    # w2(i) = (case[:w][i] / VCCF)^2 * VCCF

    # case[:KE] = [sum(u2(i)) + sum(v2(i)) + sum(w2(i)) for i in 1:Nt] ./ sum(VCCC)

    u2X(i) = (case[:uX][i] / VFCC)^2 * VFCC
    v2X(i) = (case[:vX][i] / VCFC)^2 * VCFC
    w2X(i) = (case[:wX][i] / VCCF)^2 * VCCF

    case[:KEX] = [sum(u2X(i)) + sum(v2X(i)) + sum(w2X(i)) for i in 1:Nt] ./ sum(VCCC)

    for i in 1:Nt
        case[:bX][i] ./= VCCC
    end

    return case
end

folder = "idealized_coast/"

ah = load_case(folder, "CATKE", "split_free_surface_", "QuasiAdamsBashforth2")
rh = load_case(folder, "CATKE", "split_free_surface_", "SplitRungeKutta3")
@time al = load_case(folder, "CATKE", "split_free_surface_lowres_", "QuasiAdamsBashforth2")
@time rl = load_case(folder, "CATKE", "split_free_surface_lowres_", "SplitRungeKutta3")
@time ri = load_case(folder, "CATKE", "implicit_free_surface_lowres_", "SplitRungeKutta3")


iter = Observable(1)

# xah, yah, zah = nodes(ah[:TX])
# xrh, yrh, zrh = nodes(rh[:TX])
xal, yal, zal = nodes(ri[:TX])
xrl, yrl, zrl = nodes(rl[:TX])

# bah = @lift(interior(rh[:bX][$iter], 1, :, :))
# brh = @lift(interior(ah[:bX][$iter], 1, :, :))
bal = @lift(interior(al[:bX][$iter], 1, :, :))
bri = @lift(interior(ri[:bX][$iter], 1, :, :))
brl = @lift(interior(rl[:bX][$iter], 1, :, :))

fig = Figure(size = (1000, 600))

ax = Axis(fig[1, 1]) 
# heatmap!(ax, ya, za,  u_a,  colorrange=(-60000, 60000))
# contour!(ax, yah, zah, bah, levels=range(-0.24, -0.1748, length=20), color=:blue, linewidth=2, linestyle=:dash)
# contour!(ax, yrh, zrh, brh, levels=range(-0.24, -0.1748, length=20), color=:red,  linewidth=2, linestyle=:dash)
contour!(ax, yal, zal, bal, levels=range(-0.24, -0.1748, length=20), color=:green, linewidth=2, linestyle=:solid)
contour!(ax, yal, zal, bri, levels=range(-0.24, -0.1748, length=20), color=:blue, linewidth=2, linestyle=:solid)
contour!(ax, yrl, zrl, brl, levels=range(-0.24, -0.1748, length=20), color=:red,  linewidth=2, linestyle=:solid)


add_surface_dynamics!(ri, folder, "CATKE", "implicit_free_surface_lowres_", "SplitRungeKutta3")
add_surface_dynamics!(rl, folder, "CATKE", "split_free_surface_lowres_", "SplitRungeKutta3")    
add_surface_dynamics!(al, folder, "CATKE", "split_free_surface_lowres_", "QuasiAdamsBashforth2")


# uatmp = XFaceField(usa.grid; indices = (:, :, usa.grid.Nz))
# vatmp = YFaceField(usa.grid; indices = (:, :, usa.grid.Nz))
# urtmp = XFaceField(usr.grid; indices = (:, :, usr.grid.Nz))
# vrtmp = YFaceField(usr.grid; indices = (:, :, usr.grid.Nz))

# kea = Field(@at((Center, Center, Center), uatmp^2 + vatmp^2))
# ker = Field(@at((Center, Center, Center), urtmp^2 + vrtmp^2))
# ζa  = Field(KernelFunctionOperation{Face, Face, Center}(Oceananigans.Operators.ζ₃ᶠᶠᶜ, usa.grid, uatmp, vatmp); indices = (:, :, usa.grid.Nz))
# ζr  = Field(KernelFunctionOperation{Face, Face, Center}(Oceananigans.Operators.ζ₃ᶠᶠᶜ, usr.grid, urtmp, vrtmp); indices = (:, :, usr.grid.Nz))

# Va = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Center), Oceananigans.Operators.volume, usa.grid)
# Vr = Oceananigans.AbstractOperations.grid_metric_operation((Center, Center, Center), Oceananigans.Operators.volume, usr.grid)

# us_a = @lift(interior(Field(usa[$iter] / Va), :, :, 1))
# vs_a = @lift(interior(Field(vsa[$iter] / Va), :, :, 1))
# us_r = @lift(interior(Field(usr[$iter] / Vr), :, :, 1))
# vs_r = @lift(interior(Field(vsr[$iter] / Vr), :, :, 1))
# ke_a = @lift begin
#     set!(uatmp, usa[$iter]/Va)
#     set!(vatmp, vsa[$iter]/Va)
#     Oceananigans.BoundaryConditions.fill_halo_regions!((uatmp, vatmp))
#     compute!(kea)
#     interior(kea, :, :, 1)
# end
# ke_r = @lift begin
#     set!(urtmp, usr[$iter]/Vr)
#     set!(vrtmp, vsr[$iter]/Vr)
#     Oceananigans.BoundaryConditions.fill_halo_regions!((urtmp, vrtmp))
#     compute!(ker)
#     interior(ker, :, :, 1)
# end

# ζ_a = @lift begin
#     set!(uatmp, usa[$iter]/Va)
#     set!(vatmp, vsa[$iter]/Va)
#     Oceananigans.BoundaryConditions.fill_halo_regions!((uatmp, vatmp))
#     compute!(ζa)
#     interior(ζa, :, :, 1)
# end

# ζ_r = @lift begin
#     set!(urtmp, usr[$iter]/Vr)
#     set!(vrtmp, vsr[$iter]/Vr)
#     Oceananigans.BoundaryConditions.fill_halo_regions!((urtmp, vrtmp))
#     compute!(ζr)
#     interior(ζr, :, :, 1)
# end

# xa, ya, za = nodes(usa)
# xr, yr, zr = nodes(usr)

# fig  = Figure(size = (800, 400))
# axua  = Axis(fig[1, 1])
# axva  = Axis(fig[1, 2])
# axka  = Axis(fig[2, 1])
# axza  = Axis(fig[2, 2])

# axur  = Axis(fig[3, 1])
# axvr  = Axis(fig[3, 2])
# axkr  = Axis(fig[4, 1])
# axzr  = Axis(fig[4, 2])

# heatmap!(axua, xa, ya, us_a, colorrange=(-0.5, 0.5))
# heatmap!(axva, xa, ya, vs_a, colorrange=(-0.5, 0.5))
# heatmap!(axka, xa, ya, ke_a, colorrange=(0, 0.5), colormap = :ice)
# heatmap!(axza, xa, ya, ζ_a, colorrange=(-1e-4, 1e-4), colormap = :balance)

# heatmap!(axur, xr, yr, us_r, colorrange=(-0.5, 0.5))
# heatmap!(axvr, xr, yr, vs_r, colorrange=(-0.5, 0.5))
# heatmap!(axkr, xr, yr, ke_r, colorrange=(0, 0.5), colormap = :ice)
# heatmap!(axzr, xr, yr, ζ_r, colorrange=(-1e-4, 1e-4), colormap = :balance)

# for ax in (axua, axva, axka, axza, axur, axvr, axkr, axzr)
#     hidedecorations!(ax)
# end