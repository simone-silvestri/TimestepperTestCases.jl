using Statistics: mean

function get_dissipation(filename, var)
    Ax = FieldTimeSeries(filename, "A" * var * "x")
    Ay = FieldTimeSeries(filename, "A" * var * "y")
    Az = FieldTimeSeries(filename, "A" * var * "z")

    Dz = FieldTimeSeries(filename, "D" * var * "z")

    b = FieldTimeSeries(filename, var)
    Gbx = FieldTimeSeries{Nothing, Center, Center}(grid, times)
    Gby = FieldTimeSeries{Nothing, Face,   Center}(grid, times)
    Gbz = FieldTimeSeries{Nothing, Center, Face}(grid, times)

    grid  = Abxt.grid
    times = Abxt.times
    Nx, Ny, Nz = size(grid)    
    VFCC = KernelFunctionOperation{Face, Center, Center}(Oceananigans.Operators.Vᶠᶜᶜ, grid)
    VCFC = KernelFunctionOperation{Center, Face, Center}(Oceananigans.Operators.Vᶜᶠᶜ, grid)
    VCCF = KernelFunctionOperation{Center, Center, Face}(Oceananigans.Operators.Vᶜᶜᶠ, grid)

    VFCC = compute!(Field(VFCC))
    VCFC = compute!(Field(VCFC))
    VCCF = compute!(Field(VCCF))
    
    for t in eachindex(times)
        set!(Gbx[t], mean(Gbxt[t] * VFCC, dims=1))
        set!(Gby[t], mean(Gbyt[t] * VCFC, dims=1))
        interior(Gby[t], :, 1, :) .= interior(Gby[t], :, 2, :) 
        interior(Gby[t], :, Ny+1, :) .= interior(Gby[t], :, Ny, :) 
        set!(Gbz[t], mean(Gbzt[t] * VCCF, dims=1))
        set!(Dbz[t], mean(Dbzt[t], dims=1))
    end
        
    return (; Abx, Aby, Abz, Gbx, Gby, Gbz, Dbz)
end