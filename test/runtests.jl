using TimestepperTestCases
using Test

@testset "Baroclinic stability limits" begin
    @test stability_limit((3, 2, 1)) ≈ sqrt(3) atol=1e-6
    @test stability_limit((4, 3, 2, 1)) ≈ 2sqrt(2) atol=1e-6
end

using Oceananigans

@testset "Channel energy and density diagnostics" begin
    grid = RectilinearGrid(size=(4, 4, 4), extent=(4, 4, 4), topology=(Periodic, Periodic, Periodic))
    u, v, w = XFaceField(grid), YFaceField(grid), ZFaceField(grid)
    volumes = (CenterField(grid), XFaceField(grid), YFaceField(grid), ZFaceField(grid))
    set!(u, 2)
    set!(v, 3)
    set!(w, 6)
    foreach(volume -> set!(volume, 1), volumes)
    case = Dict(:u => [u], :v => [v], :w => [w],
                :VCCC => [volumes[1]], :VFCC => [volumes[2]], :VCFC => [volumes[3]], :VCCF => [volumes[4]])
    @test TimestepperTestCases.KineticEnergy(case, 1) ≈ (2^2 + 3^2 + 6^2) / 2
    @test TimestepperTestCases.MeanKineticEnergy(case, 1) ≈ (2^2 + 3^2 + 6^2) / 2
    b = CenterField(grid)
    set!(b, 0.01)
    density = TimestepperTestCases.DensityField(b; ρ₀=1025, g=10)
    @test all(interior(density) .≈ 1025 * (1 - 0.01 / 10))
end
