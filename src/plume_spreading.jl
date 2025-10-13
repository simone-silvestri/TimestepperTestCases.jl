using Oceananigans

plume_spreading_timestep(::Val{:QuasiAdamsBashforth2}) = 20minutes
plume_spreading_timestep(::Val{:SplitRungeKutta3})     = 60minutes
plume_spreading_timestep(::Val{:SplitRungeKutta6})     = 120minutes

function plume_spreading(timestepper::Symbol)

    grid = RectilinearGrid(size = (48, 48, 8),
                           x = (-500kilometers, 500kilometers),
                           y = (-500kilometers, 500kilometers),
                           z = (-1kilometers, 0),
                           halo = (6, 6, 6),
                           topology = (Periodic, Periodic, Bounded))


end