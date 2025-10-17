using Statistics: mean

function get_dissipation(filename, var)
    Ax = FieldTimeSeries(filename, "A" * var * "x")
    Az = FieldTimeSeries(filename, "A" * var * "z")
    Gx = FieldTimeSeries(filename, "G" * var * "x")
    Gz = FieldTimeSeries(filename, "G" * var * "z")
    Vx = FieldTimeSeries(filename, "VFC")
    Vt = FieldTimeSeries(filename, "VCC")
    Vz = FieldTimeSeries(filename, "VCF")

    grid  = Ax.grid
    times = Ax.times
    Nx, Ny, Nz = size(grid)    
    
    for t in eachindex(times)
        set!(Gx[t], Gx[t] * Vx)
        set!(Gz[t], Gz[t] * Vz)
    end

    return (; Ax, Az, Gx, Gz)
end