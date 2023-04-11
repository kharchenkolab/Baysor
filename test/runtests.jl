using Test
using DataFrames
using Distributions
using LinearAlgebra
using Statistics
using StatsBase

import Random: seed!
import Baysor: BPR, DAT, REP, Utils

module DataWrappers

    using DataFrames
    import Baysor: BPR

    function get_bmm_data(; n_mols::Int=5000, confidence::Bool=false, scale::Float64=0.1, min_molecules_per_cell::Int=10, do_maximize::Bool=false, kwargs...)
        df = DataFrame(
            :x => rand(n_mols), :y => rand(n_mols),
            :gene => rand(1:10, n_mols),
            :confidence => confidence ? ones(n_mols) : rand(n_mols)
        )
        adj_list = BPR.build_molecule_graph(
            df; use_local_gene_similarities=false, adjacency_type=:triangulation
        )
        bmd = BPR.initialize_bmm_data(
            df; scale, min_molecules_per_cell, n_cells_init=(n_mols ÷ 5), adj_list, kwargs...
        );
        if do_maximize
            BPR.maximize!(bmd)
        end

        return bmd
    end

end

@testset "Baysor" begin
    @testset "utils" begin
        @testset "convex_hull" begin
            chull = BPR.convex_hull([[0, 0], [1, 0], [0, 0], [1, 1], [0, 1]])
            @test all(chull .== [0 0 1 1 0; 0 1 1 0 0])
            @test BPR.area(chull) ≈ 1.0

            chull = BPR.convex_hull([[0, 0], [1, 0], [0, 0], [1, 1], [0, 1], [0, 2], [2, 0]])
            @test all(chull .== [0 0 2 0; 0 2 0 0])
            @test BPR.area(chull) ≈ 2.0

            chull = BPR.convex_hull([[0, 0], [1, 0], [0, 0], [1, 1], [0, 1], [0, 2], [2, 0], [2, 2]])
            @test all(chull .== [0 0 2 2 0; 0 2 2 0 0])
            @test BPR.area(chull) ≈ 4.0
        end

        @testset "encode_genes" begin
            genes = ["1", "2", "3"]
            gene_ids, gene_names = DAT.encode_genes(genes)
            @test all(gene_ids .== [1, 2, 3])
            @test all(gene_names .== genes)

            genes = ["1", "2", missing, "3", missing]
            gene_ids, gene_names = DAT.encode_genes(genes)
            @test all(skipmissing(gene_ids) .== [1, 2, 3])
            @test ismissing(gene_ids[3]) && ismissing(gene_ids[5])
        end

        @testset "other" begin
            @test Utils.get_cell_name(0) == ""
            @test Utils.get_cell_name(1; run_id="") == "1"
            @test Utils.get_cell_name(2; run_id="RI") == "CRI-2"

            @test Utils.isnoise(Utils.get_cell_name(0))
            @test Utils.isnoise(Utils.get_cell_name(0; run_id="RI"))

            @test all(Utils.isnoise.([0, 1, 2, 1]) .== [true, false, false, false])
            @test all(Utils.isnoise.(["", "0", "2", "1"]) .== [true, false, false, false])
        end
    end

    @testset "distributions" begin
        @testset "ShapePrior" begin
            seed!(42)

            means = BPR.MeanVec{2}([10.0, 20.0])
            std_stds = BPR.MeanVec{2}([5.0, 10.0])
            stds = hcat([Vector(BPR.sample_var(BPR.ShapePrior{2}(means, std_stds, 1000))) for i in 1:100000]...) .^ 0.5
            @test all(stds .> 0)
            @test all(abs.(vec(mean(stds, dims=2)) .- means) .< 1)

            means = BPR.MeanVec{2}([1000.0, 2000.0])
            stds = hcat([Vector(BPR.sample_var(BPR.ShapePrior{2}(means, std_stds, 1000))) for i in 1:100000]...) .^ 0.5
            @test all(stds .> 0)
            @test all(abs.(vec(mean(stds, dims=2)) .- means) .< 1)
            @test all(abs.(vec(std(stds, dims=2)) .- std_stds) .< 1)
        end

        @testset "CategoricalSmoothed" begin
            dist = BPR.CategoricalSmoothed([0.0, 0.0, 0.0])
            BPR.maximize!(dist, [missing, 1, 2], [0.5, 0.5, 0.5])
            @test all(dist.counts .≈ [0.5, 0.5, 0.0])
            @test dist.sum_counts ≈ 1.0
        end
    end

    @testset "data_processing" begin
        @testset "initialization" begin
            n_mols = 1000
            df = DataFrame(:x => rand(n_mols), :y => rand(n_mols), :gene => rand(1:10, n_mols))

            df[!, :cluster] = rand(1:5, n_mols)
            df[!, :prior_segmentation] = rand(1:100, n_mols)

            adj_list = BPR.build_molecule_graph(
                df; use_local_gene_similarities=false, adjacency_type=:triangulation
            )

            bm_data = BPR.initialize_bmm_data(
                df; scale=6.0, min_molecules_per_cell=30, n_cells_init=200, adj_list
            );
            @test all(bm_data.cluster_per_molecule .== df.cluster)
            @test all(bm_data.segment_per_molecule .== df.prior_segmentation)

            @test all(length(bm_data.adj_list.ids[i]) == length(bm_data.adj_list.weights[i]) for i in eachindex(bm_data.adj_list))
            @test !any(i in bm_data.adj_list.ids[i] for i in eachindex(bm_data.adj_list))

            for i in 5:10:55
                init_params = BPR.cell_centers_uniformly(df, i; scale=10.0)
                @test size(init_params.centers, 1) == i
                @test size(init_params.centers, 2) == 2
                @test length(init_params.covs) == i
                @test maximum(init_params.assignment) == i
            end
        end

        @testset "boundary" begin
            @test all(BPR.border_edges_to_poly([[1, 1]]) .== [1, 1])
            @test all(BPR.border_edges_to_poly([[1, 2], [2, 1]]) .== [2, 1, 2])
            @test all(BPR.border_edges_to_poly([[1, 2], [3, 1], [2, 3]]) .== [2, 3, 1, 2])
        end

        @testset "parse_parameters" begin
            @test BPR.parse_scale_std("23.55%", 121.0) ≈ 121.0 * 0.2355
            @test BPR.parse_scale_std(33.87, 121.0) ≈ 33.87
            @test BPR.parse_scale_std(nothing, 121.0) ≈ 121.0 * 0.25
        end

        @testset "molecule_graph" begin
            for adj_type in [:triangulation, :knn, :both]
                adj_list = BPR.build_molecule_graph(DataFrame(rand(1000, 2), [:x, :y]); adjacency_type=adj_type, k_adj=30);
                @test all(length.(adj_list.ids) .== length.(adj_list.weights))
                @test all(length.(adj_list.ids) .== length.(unique.(adj_list.ids)))
            end
        end

        @testset "filter_segmentation" begin
            labels = [1 2 3; 1 2 4; 5 2 7];
            position_data = hcat(vec([[id[1], id[2]] for id in CartesianIndices(labels)])...)
            df_spatial = DataFrame(:x => position_data[1,:], :y => position_data[2,:])
            assignment = DAT.staining_value_per_transcript(df_spatial, labels)

            expected_results = [
                labels,
                labels,
                [1 2 0; 1 2 0; 0 2 0],
                [0 2 0; 0 2 0; 0 2 0],
                [0 0 0; 0 0 0; 0 0 0]
            ]
            for m in 0:4
                l, ta = DAT.filter_segmentation_labels!(deepcopy(labels), deepcopy(assignment), min_molecules_per_segment=m)
                @test all(l .== expected_results[m + 1])

                a = DAT.filter_segmentation_labels!(deepcopy(assignment), min_molecules_per_segment=m)
                @test all(a .== DAT.staining_value_per_transcript(df_spatial, l))
                @test all(a .== ta)
            end
        end
    end

    @testset "bmm_algorithm" begin
        @testset "noise_composition_density" begin
            for conf in [false, true]
                bm_data = DataWrappers.get_bmm_data(confidence=conf, do_maximize=true);
                dens_exp = mean([mean(c.composition_params.counts[c.composition_params.counts .> 0] ./ c.composition_params.sum_counts) for c in bm_data.components]);
                dens_obs = BPR.noise_composition_density(bm_data)
                @test abs(dens_exp - dens_obs) < 1e-10
            end
        end

        @testset "distribution_sampling" begin
            bm_data = DataWrappers.get_bmm_data(confidence=true, do_maximize=true);
            @test_nowarn for i in 1:1000
                BPR.sample_distribution!(bm_data)
            end
        end

        @testset "drop_unused_components" begin
            bm_data = DataWrappers.get_bmm_data(n_mols=1000, confidence=true);
            bm_data.assignment[rand(1:length(bm_data.assignment), 100)] .= 0
            BPR.maximize!(bm_data)
            BPR.drop_unused_components!(bm_data)

            @test all(BPR.num_of_molecules_per_cell(bm_data) .== [c.n_samples for c in bm_data.components])
        end

        @testset "maximize_mvnormal" begin
            is_not_nan = map(1:10000) do i
                n_samps = rand(0:100)
                d = BPR.MvNormalF([0., 0.], diagm(0 => ones(2)));
                BPR.maximize!(d, rand(2, n_samps); center_probs=rand(n_samps) .* rand(Binomial(1, 0.5), n_samps))
                !isnan(d.Σ[1])
            end

            @test all(is_not_nan)
        end

        @testset "synthetic_run" begin
            seed!(42)
            n_components = 10;
            n_genes = 20
            scale = 0.2;
            frame_size = (100, 200)
            cell_size = 50
            noise_size = 100
            report_file = tempname()

            for it in 1:6
                centers = hcat(rand(n_components) * frame_size[1], rand(n_components) * frame_size[2]);
                sizes = rand(Poisson(cell_size), n_components);
                expressions = rand(Dirichlet(n_genes, 0.1), n_components);

                df_spatial = vcat([DataFrame(
                    :x => rand(Normal(centers[i,1], scale), sizes[i]),
                    :y => rand(Normal(centers[i,2], scale), sizes[i]),
                    :gene => rand(Categorical(expressions[:,i]), sizes[i]),
                    :cell => i
                ) for i in 1:n_components]...)

                df_spatial = vcat(df_spatial, DataFrame(
                    :x => rand(noise_size) * frame_size[1],
                    :y => rand(noise_size) * frame_size[2],
                    :gene => rand(1:n_genes, noise_size),
                    :cell => 0
                ));

                if it % 2 == 0
                    df_spatial[!, :z] = rand(1:3, size(df_spatial, 1)) .* scale ./ 10
                end

                if it > 4
                    df_spatial[!, :cluster] = df_spatial.cell .+ 1
                end

                adj_list = BPR.build_molecule_graph(df_spatial)
                bm_data = BPR.initialize_bmm_data(
                    df_spatial; scale_std="5%", min_molecules_per_cell=10, n_cells_init=size(df_spatial, 1) ÷ 3,
                    adj_list, scale
                );
                BPR.bmm!(bm_data; n_iters=350, new_component_frac=0.3, min_molecules_per_cell=10, assignment_history_depth=100, verbose=false);

                conj_table = counts(BPR.estimate_assignment_by_history(bm_data)[1], df_spatial.cell)
                @test all([all(vec(mapslices(maximum, conj_table ./ sum(conj_table, dims=d), dims=d)) .> 0.8) for d in 1:2])
                @test maximum(conj_table[2:end,1]) <= 3 # If there are new cells consisting of noise, they should be very small

                # Test that it works
                seg_df = BPR.get_segmentation_df(bm_data)
                gene_names = ["g$i" for i in 1:maximum(df_spatial.gene)]
                seg_df = BPR.get_segmentation_df(bm_data, gene_names)
                df_res = BPR.get_cell_stat_df(bm_data, seg_df)
                cm1 = BPR.convert_segmentation_to_counts(BPR.composition_data(bm_data), bm_data.assignment)
                cm2 = BPR.convert_segmentation_to_counts(BPR.composition_data(bm_data), bm_data.assignment, gene_names)
                @test all((size(cm2[:, 2:end]) .== size(cm1)))

                segmented_df, cell_stat_df, cm = BPR.get_segmentation_results(bm_data, gene_names; run_id="test")
                @test length(unique(cell_stat_df.cell)) == size(cm, 2)
                @test length(cell_stat_df.cell) .== length(unique(cell_stat_df.cell))

                cs1,cs2 = (sort(unique(vs)) for vs in (segmented_df.cell, cell_stat_df.cell));
                cs1 = cs1[cs1 .!= ""]

                @test length(cs1) == length(cs2)
                @test all(cs1 .== cs2)

                # Test that reports work
                REP.plot_diagnostics_panel(segmented_df, segmented_df.cell, bm_data.tracer; file=report_file)
            end
        end
    end

    @testset "reports" begin
        @testset "plot_molecules" begin
            # Just test that it works. Also needed for pre-compilation derectives
            n_obs = 100
            n_cells = 10
            p_df = DataFrame(
                :x => (1:n_obs) ./ 10, :y => rand(n_obs),
                :gene => rand(1:10, n_obs),
                :cell => repeat(1:n_cells, inner=(n_obs ÷ n_cells))
            );

            polygons = BPR.boundary_polygons(BPR.position_data(p_df), p_df.cell)
            @test length(polygons) == 10

            fig = REP.plot_molecules(p_df);
            fig = REP.plot_molecules(p_df, polygons, annotation=:cell);
            REP.makie_to_base64(fig);
        end
    end
end