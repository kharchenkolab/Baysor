using Baysor, Test

@testset "Baysor" begin
    @testset "utils" begin
        @testset "convex_hull" begin
            chull = Baysor.convex_hull([[0, 0], [1, 0], [0, 0], [1, 1], [0, 1]])
            @test all(chull .== [0 0 1 1 0; 0 1 1 0 0])
            @test Baysor.area(chull) ≈ 1.0

            chull = Baysor.convex_hull([[0, 0], [1, 0], [0, 0], [1, 1], [0, 1], [0, 2], [2, 0]])
            @test all(chull .== [0 0 2 0; 0 2 0 0])
            @test Baysor.area(chull) ≈ 2.0

            chull = Baysor.convex_hull([[0, 0], [1, 0], [0, 0], [1, 1], [0, 1], [0, 2], [2, 0], [2, 2]])
            @test all(chull .== [0 0 2 2 0; 0 2 2 0 0])
            @test Baysor.area(chull) ≈ 4.0
        end

        @testset "utils" begin
            @test all([Baysor.interpolate_linear(x, 0.0, 1.0; y_start=0.0, y_end=1.0) ≈ x for x in range(0, 1.0, length=50)])
        end
    end
end