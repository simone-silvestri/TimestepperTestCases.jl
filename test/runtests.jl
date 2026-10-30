using TimestepperTestCases
using Test

@testset "Baroclinic stability limits" begin
    @test stability_limit((3, 2, 1)) ≈ sqrt(3) atol=1e-6
    @test stability_limit((4, 3, 2, 1)) ≈ 2sqrt(2) atol=1e-6
end
