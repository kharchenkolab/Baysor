using Baysor, Test
using DataFrames
using Statistics
import Random: seed!

B = Baysor

@testset "Baysor" begin
    @testset "utils" begin
        @testset "convex_hull" begin
            chull = B.convex_hull([[0, 0], [1, 0], [0, 0], [1, 1], [0, 1]])
            @test all(chull .== [0 0 1 1 0; 0 1 1 0 0])
            @test B.area(chull) ≈ 1.0

            chull = B.convex_hull([[0, 0], [1, 0], [0, 0], [1, 1], [0, 1], [0, 2], [2, 0]])
            @test all(chull .== [0 0 2 0; 0 2 0 0])
            @test B.area(chull) ≈ 2.0

            chull = B.convex_hull([[0, 0], [1, 0], [0, 0], [1, 1], [0, 1], [0, 2], [2, 0], [2, 2]])
            @test all(chull .== [0 0 2 2 0; 0 2 2 0 0])
            @test B.area(chull) ≈ 4.0
        end

        @testset "utils" begin
            @test all([B.interpolate_linear(x, 0.0, 1.0; y_start=0.0, y_end=1.0) ≈ x for x in range(0, 1.0, length=50)])
        end
    end

    @testset "distributions" begin
        @testset "ShapePrior" begin
            seed!(42)

            means = B.MeanVec([10.0, 20.0])
            std_stds = B.MeanVec([5.0, 10.0])
            stds = hcat([Vector(B.sample_var(B.ShapePrior(means, std_stds, 1000))) for i in 1:100000]...) .^ 0.5
            @test all(stds .> 0)
            @test all(abs.(vec(mean(stds, dims=2)) .- means) .< 1)

            means = B.MeanVec([1000.0, 2000.0])
            stds = hcat([Vector(B.sample_var(B.ShapePrior(means, std_stds, 1000))) for i in 1:100000]...) .^ 0.5
            @test all(stds .> 0)
            @test all(abs.(vec(mean(stds, dims=2)) .- means) .< 1)
            @test all(abs.(vec(std(stds, dims=2)) .- std_stds) .< 1)
        end
    end

    @testset "data_processing" begin
        @testset "initialization" begin
            df = DataFrame(:x => rand(1000), :y => rand(1000), :gene => rand(1:10, 1000))

            for i in 1:10
                bm_data_arr = B.initial_distribution_arr(df, n_frames=i, scale=6.0, min_molecules_per_cell=30);
                @test length(bm_data_arr) <= i
            end

            df[!, :cluster] = rand(1:5, 1000)
            df[!, :prior_segmentation] = rand(1:100, 1000)
            bm_data = B.initial_distribution_arr(df, n_frames=1, scale=6.0, min_molecules_per_cell=30)[1];
            @test all(bm_data.cluster_per_molecule .== df.cluster)
            @test all(bm_data.segment_per_molecule .== df.prior_segmentation)
        end

        @testset "parse_parameters" begin
            @test B.parse_scale_std("23.55%", 121.0) ≈ 121.0 * 0.2355
            @test B.parse_scale_std(33.87, 121.0) ≈ 33.87
            @test B.parse_scale_std(nothing, 121.0) ≈ 121.0 * 0.25
        end

        @testset "molecule_graph" begin
            for adj_type in [:triangulation, :knn, :both]
                adj_points, adj_weights = B.build_molecule_graph(DataFrame(rand(1000, 2), [:x, :y]); adjacency_type=adj_type, k_adj=30)[1:2];
                @test all(length.(adj_points) .== length.(adj_weights))
                @test all(length.(adj_points) .== length.(unique.(adj_points)))
            end
        end
    end
end