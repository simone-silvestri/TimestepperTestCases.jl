using Oceananigans.BoundaryConditions
using Oceananigans.Operators

u2(case, i) = (case[:u][i])^2 / case[:VFCC][i]
v2(case, i) = (case[:v][i])^2 / case[:VCFC][i]
w2(case, i) = (case[:w][i])^2 / case[:VCCF][i]

um2(case, i) = (mean(case[:u][i], dims=1))^2 / mean(case[:VFCC][i], dims=1)
vm2(case, i) = (mean(case[:v][i], dims=1))^2 / mean(case[:VCFC][i], dims=1)
wm2(case, i) = (mean(case[:w][i], dims=1))^2 / mean(case[:VCCF][i], dims=1)

Sm(case, i)   = mean(case[:S][i] / case[:VCCC][i], dims=1)

function save_variable(name, location, path, grid, times)
    new_path = path * "_full.jld2"
    old_path = path * ".jld2"
    ftsn = FieldTimeSeries{location...}(grid, times; backend = OnDisk(), path=new_path, name)
    ftso = FieldTimeSeries(old_path, name)

    for t in 1:length(times)
        @info "Saving $name time step $t / $(length(times))"
        set!(ftsn, ftso[t], t)
    end
end

const FCC = (Face, Center, Center)
const CFC = (Center, Face, Center)
const CCF = (Center, Center, Face)
const CCC = (Center, Center, Center)
const CCN = (Center, Center, Nothing)

function save_case(folder, closure, suffix, timestepper)
    path  = folder * "idealized_coast_" * suffix * timestepper * "_" * closure 
    grid  = FieldTimeSeries(path * ".jld2", "u").grid
    times = FieldTimeSeries(path * ".jld2", "u").times
    args  = (path, grid, times)
    save_variable("u", FCC, args...)
    save_variable("v", CFC, args...)
    save_variable("w", CCF, args...)
    save_variable("S", CCC, args...)
    save_variable("T", CCC, args...)
    save_variable("b", CCC, args...)
    save_variable("η", CCN, args...)
    
    save_variable("Abx", FCC, args...)
    save_variable("Aby", CFC, args...)
    save_variable("Abz", CCF, args...)
    save_variable("Gbx", FCC, args...)
    save_variable("Gby", CFC, args...)
    save_variable("Gbz", CCF, args...)

    save_variable("ASx", FCC, args...)
    save_variable("ASy", CFC, args...)
    save_variable("ASz", CCF, args...)
    save_variable("GSx", FCC, args...)
    save_variable("GSy", CFC, args...)
    save_variable("GSz", CCF, args...)

    save_variable("ATx", FCC, args...)
    save_variable("ATy", CFC, args...)
    save_variable("ATz", CCF, args...)
    save_variable("GTx", FCC, args...)
    save_variable("GTy", CFC, args...)
    save_variable("GTz", CCF, args...)
end

function load_idealized_coast(folder, closure, suffix, timestepper)
    path = folder * "idealized_coast_" * suffix * timestepper * "_" * closure * ".jld2"
    case = Dict()

    case[:u] = FieldTimeSeries(path, "u"; backend=OnDisk())
    case[:v] = FieldTimeSeries(path, "v"; backend=OnDisk())
    case[:w] = FieldTimeSeries(path, "w"; backend=OnDisk())
    case[:S] = FieldTimeSeries(path, "S"; backend=OnDisk())
    case[:T] = FieldTimeSeries(path, "T"; backend=OnDisk())
    case[:b] = FieldTimeSeries(path, "b"; backend=OnDisk())
    case[:η] = FieldTimeSeries(path, "η")

    fill_halo_regions!(case[:η])

    grid = case[:u].grid
    times = case[:u].times

    Nx, Ny, Nz = size(grid)

    VCCC = FieldTimeSeries{Center, Center, Center}(grid, times)
    VFCC = FieldTimeSeries{Face,   Center, Center}(grid, times)
    VCFC = FieldTimeSeries{Center, Face,   Center}(grid, times)
    VCCF = FieldTimeSeries{Center, Center, Face  }(grid, times)

    Nt = length(times)

    params = Oceananigans.Utils.KernelParameters(0:Nx+1, 0:Ny+1, 0:Nz+1)
    _compute_volumes_kernel! = Oceananigans.Utils.configure_kernel(CPU(), grid, params, _compute_volumes!)[1]

    for t in 1:Nt
        _compute_volumes_kernel!(VCCC[t], VFCC[t], VCFC[t], VCCF[t], grid, case[:η][t])
    end
    
    fill_halo_regions!(VCCC)
    fill_halo_regions!(VFCC)
    fill_halo_regions!(VCFC)
    fill_halo_regions!(VCCF)

    case[:VCCC] = VCCC
    case[:VFCC] = VFCC
    case[:VCFC] = VCFC
    case[:VCCF] = VCCF
    GC.gc()

    case[:Abx] = FieldTimeSeries(path, "Abx"; backend=OnDisk())
    case[:Aby] = FieldTimeSeries(path, "Aby"; backend=OnDisk())
    case[:Abz] = FieldTimeSeries(path, "Abz"; backend=OnDisk())

    case[:abx] = [sum(case[:Abx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:aby] = [sum(case[:Aby][i]) for i in 1:Nt] ./ [sum(case[:VCFC][i]) for i in 1:Nt]
    case[:abz] = [sum(case[:Abz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]
    case[:abt] = case[:abx] .+ case[:aby] .+ case[:abz]

    case[:Gbx] = FieldTimeSeries(path, "Gbx"; backend=OnDisk())
    case[:Gby] = FieldTimeSeries(path, "Gby"; backend=OnDisk())
    case[:Gbz] = FieldTimeSeries(path, "Gbz"; backend=OnDisk())

    case[:gbx] = [sum(case[:Gbx][i]) for i in 1:Nt] ./ [sum(case[:VFCC][i]) for i in 1:Nt]
    case[:gby] = [sum(case[:Gby][i]) for i in 1:Nt] ./ [sum(case[:VCFC][i]) for i in 1:Nt]
    case[:gbz] = [sum(case[:Gbz][i]) for i in 1:Nt] ./ [sum(case[:VCCF][i]) for i in 1:Nt]    
    case[:gbt] = case[:gbx] .+ case[:gby] .+ case[:gbz]
    GC.gc()

    case[:KE]  = [sum(u2(case, i))  + sum(v2(case, i))  + sum(w2(case, i))  for i in 1:Nt] ./ [sum(case[:VCCC][i]) for i in 1:Nt]
    case[:MKE] = [sum(um2(case, i)) + sum(vm2(case, i)) + sum(wm2(case, i)) for i in 1:Nt] ./ [sum(mean(case[:VCCC][i], dims=1)) for i in 1:Nt]
    case[:η2]  = [mean(case[:η][i]^2) for i in 1:Nt]

    GC.gc()

    return case
end

using FFTW

function power_spectrum_x(var, x)
    Nx      = length(x)
    Nfx     = Int64(Nx)
    spectra = zeros(ComplexF64, Int(Nfx/2))
    dx      = x[2] - x[1]

    freqs = fftfreq(Nfx, 1.0 / dx) # 0,+ve freq,-ve freqs (lowest to highest)
    freqs = freqs[1:Int(Nfx/2)] .* 2.0 .* π

    fourier      = fft(var) / Nfx
    spectra[1]  += fourier[1] .* conj(fourier[1])

    for m in 2:Int(Nfx/2)
        spectra[m] += 2.0 * fourier[m] * conj(fourier[m]) # factor 2 for neg freq contribution
    end

    return spectra, freqs
end

function y_average_spectra(var::Field, irange, jrange; k=40, spectra=power_spectrum_x)

    xdomain = xnodes(var)[irange]
    ydomain = ynodes(var)[jrange]
    spec, freq = spectra(interior(var, irange, jrange[1], k), xdomain)

    for j in jrange[2:end]
        spec .+= spectra(interior(var, irange, j, k), xdomain)[1]
    end

    spec .*= 1 / length(jrange)

    return spec, freq
end