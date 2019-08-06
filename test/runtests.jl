using Baysor, Test
using Statistics
import Random: seed!

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

    @testset "distributions" begin
        @testset "ScaledInverseChisq" begin
            seed!(42)
            for (μ, σ) in zip(rand(20) .* 1000 .+ 1000, rand(20) .* 100)
                c_sample = rand(Baysor.ScaledInverseChisq(μ=μ, σ²=σ^2), 10000000)

                @test isapprox(mean(c_sample), μ, atol=0.5)
                @test isapprox(std(c_sample), σ, atol=0.5)
            end
        end

        @testset "ShapePrior" begin
            seed!(42)

            means = [10.0, 20.0]
            std_stds = [5.0, 10.0]
            stds = hcat([Baysor.sample_var(Baysor.ShapePrior(means, std_stds)) for i in 1:100000]...) .^ 0.5
            @test all(stds .> 0)
            @test all(abs.(vec(mean(stds, dims=2)) .- means) .< 1)

            means = [1000.0, 2000.0]
            stds = hcat([Baysor.sample_var(Baysor.ShapePrior(means, std_stds)) for i in 1:100000]...) .^ 0.5
            @test all(stds .> 0)
            @test all(abs.(vec(mean(stds, dims=2)) .- means) .< 1)
            @test all(abs.(vec(std(stds, dims=2)) .- std_stds) .< 1)
        end
    end

    @testset "data_processing" begin
        @testset "initialization"
            @test parse_scale_std("23.55%", 121.0) ≈ 23.55
            @test parse_scale_std(33.87, 121.0) ≈ 33.87
            @test parse_scale_std(nothing, 121.0) ≈ 121.0 * 0.25
        end
    end
end