using Test
using DataFrames
using Distributions
using LinearAlgebra
using Statistics
using StatsBase
using StaticArrays

import Random: seed!
import Baysor: BPR, DAT, REP, Utils

module DataWrappers

    using DataFrames
    import Baysor: BPR

    function get_spatial_df(;
            n_mols::Int=5000, confidence::Bool=false, cluster::Bool=false, prior_segmentation::Bool=false, cell::Bool=false,
            n_cells=100
        )
        df = DataFrame(:x => rand(n_mols), :y => rand(n_mols), :gene => rand(1:10, n_mols))
        if confidence
            df[!, :confidence] = rand(n_mols)
        end

        if cluster
            df[!, :cluster] = rand(1:5, n_mols)
        end

        if prior_segmentation
            df[!, :prior_segmentation] = rand(1:100, n_mols)
        end

        if cell
            df[!, :cell] = rand(1:n_cells, n_mols)
        end

        return df
    end

    function get_bmm_data(; n_mols::Int=5000, confidence::Bool=false, scale::Float64=0.1, min_molecules_per_cell::Int=10, do_maximize::Bool=false, kwargs...)
        df = get_spatial_df(; n_mols, confidence=confidence)
        if !confidence
            df[!, :confidence] = ones(n_mols)
        end

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
        @testset "CategoricalSmoothed" begin
            dist = BPR.CategoricalSmoothed([0.0, 0.0, 0.0])
            BPR.maximize!(dist, [missing, 1, 2], [0.5, 0.5, 0.5])
            @test all(dist.counts .≈ [0.5, 0.5, 0.0])
            @test dist.sum_counts ≈ 1.0
            @test dist.n_genes == 2
        end
    end

    @testset "data_loading" begin
        @testset "prior" begin
            df = DataWrappers.get_spatial_df(; n_mols=1000000, n_cells=10000, confidence=false, cell=true)
            pos_data = BPR.position_data(df)
            scale, std = DAT.estimate_scale_from_assignment(pos_data, df.cell, min_molecules_per_cell=1);
            @test (scale - 0.00026) < 0.00002
            @test (std - 0.00017) < 0.00002
        end

        @testset "not_enough_molecules" begin
            df = DataWrappers.get_spatial_df(prior_segmentation=true);

            @test_throws ErrorException DAT.load_prior_segmentation!(
                ":prior_segmentation", df, min_molecules_per_segment=1000,
                min_molecules_per_cell=1000, estimate_scale=true, unassigned_label="0"
            )

            prior_seg, scale, scale_std = DAT.load_prior_segmentation!(
                ":prior_segmentation", df, min_molecules_per_segment=10,
                min_molecules_per_cell=1000, estimate_scale=false, unassigned_label="0"
            )

            @test prior_seg === nothing
            @test scale ≈ 0.0
            @test scale_std ≈ 0.0
        end
    end

    @testset "data_processing" begin
        @testset "initialization" begin
            df = DataWrappers.get_spatial_df(; n_mols=1000, confidence=false, cluster=true, prior_segmentation=true)

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
            @testset "get_n_points_in_triangle" begin
                tri = SMatrix{2, 3, Float64}([0.0 0.0; 20.0 0.0; 0.0 20.0]')

                for np in [1, 5, 10, 20]
                    cp = hcat(rand(np), rand(np))'
                    @test BPR.get_n_points_in_triangle(cp, tri) == np
                end

                @test BPR.get_n_points_in_triangle([100.; 100.;;], tri) == 0
            end

            @testset "border_edges_to_poly" begin
                @test isempty(BPR.border_edges_to_poly([(1, 1)]))
                @test isempty(BPR.border_edges_to_poly([(1, 2), (2, 1)]))
                @test all(BPR.border_edges_to_poly([(1, 2), (3, 1), (2, 3)]) .== [2, 1, 3, 2])

                cb = BPR.border_edges_to_poly([(313, 353), (353, 416), (232, 313), (232, 416)])
                @test all(cb .== [313, 353, 416, 232, 313])

                cb = BPR.border_edges_to_poly([(11, 61), (259, 264), (63, 264), (61, 63), (11, 259)])
                @test all(cb .== [61, 11, 259, 264, 63, 61])
            end

            @testset "extract_ids_per_bbox" begin
                tx, ty = [rand(10000) for _ in 1:2];
                dvals = 0.1:0.1:0.5
                bbs = [[rand() * (1 - d), rand() * (1 - d)] for d in dvals];
                bbs = [[(tx, tx + d), (ty, ty + d)] for ((tx,ty),d) in zip(bbs, dvals)];

                tids_per_bb = BPR.extract_ids_per_bbox(tx, ty, bbs)
                for (((bxs, bxe), (bys, bye)), tids) in zip(bbs, tids_per_bb)
                    @test all((tx[tids] .< bxe) .& (tx[tids] .> bxs))
                    @test all((ty[tids] .< bye) .& (ty[tids] .> bys))
                end
            end

            @testset "bounding_boxes" begin
                df_spatial = DataWrappers.get_spatial_df(; n_mols=100000, cell=true)
                pos_data = BPR.position_data(df_spatial)
                bbox_per_cell = BPR.get_boundary_box_per_cell(pos_data, df_spatial.cell)
                ids_per_bbox = BPR.extract_ids_per_bbox(df_spatial.x, df_spatial.y, bbox_per_cell)
                for (ci,ids) in enumerate(ids_per_bbox)
                    cids = findall(df_spatial.cell .== ci)
                    @test length(intersect(ids, cids)) == length(cids)
                end
            end

            @testset "boundary_polygons" begin
                df_spatial = DataWrappers.get_spatial_df(; n_mols=100000, cell=true)
                pos_data = BPR.position_data(df_spatial)
                polygons = BPR.boundary_polygons(pos_data, df_spatial.cell)

                # Basic tests
                @test length(polygons) == maximum(df_spatial.cell)
                @test all(BPR.area(Matrix(p')) > 1e-5 for (i,p) in polygons)

                # Z-stack format
                n_stacks = 3
                z_vals = rand(1:n_stacks, size(df_spatial, 1))
                pos_data = vcat(pos_data, z_vals')
                polygons = BPR.boundary_polygons_auto(pos_data, df_spatial.cell; estimate_per_z=true)[2]
                @test length(polygons) == (n_stacks + 1)
                @test "2d" in keys(polygons)
                @test all("$(Float64(v))" in keys(polygons) for v in unique(z_vals))

                n_stacks = 50
                max_slices = 10
                z_vals = rand(1:n_stacks, size(df_spatial, 1))
                pos_data[3,:] .= z_vals
                polygons = BPR.boundary_polygons_auto(pos_data, df_spatial.cell; estimate_per_z=true, max_z_slices=max_slices)[2]
                @test length(polygons) == (max_slices + 1)
                @test "2d" in keys(polygons)

                # Process 1-point cells
                n_points = 100
                pos_data = rand(2, n_points)
                polygons = BPR.boundary_polygons(pos_data, collect(1:n_points))
                @test length(polygons) == n_points
                @test all(size.(values(polygons), 1) .== 4)
            end
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

        @testset "neighborhood_composition" begin
            df = DataWrappers.get_spatial_df(; n_mols=100000, confidence=true)
            pos_data = BPR.position_data(df)
            k = 10

            for conf in [nothing, df.confidence]
                for nbd in [false, true]
                    for (norm, s) in [(true, 1.0), (false, k)]
                        cm = BPR.neighborhood_count_matrix(pos_data, df.gene, k; confidences=conf, normalize_by_dist=nbd, normalize=norm)
                        @test all(abs.(sum(cm, dims=1)[:] .- s) .< 1e-4)
                    end
                end
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

        @testset "drop_unused_components" begin
            bm_data = DataWrappers.get_bmm_data(n_mols=1000, confidence=true);
            bm_data.assignment[rand(1:length(bm_data.assignment), 100)] .= 0
            BPR.maximize!(bm_data)
            BPR.drop_unused_components!(bm_data)

            @test all(BPR.num_of_molecules_per_cell(bm_data) .== [c.n_samples for c in bm_data.components])
        end

        @testset "split_connected_components" begin
            tbd = DataWrappers.get_bmm_data(n_mols=20000);
            tbd.assignment[5001:end] .= mod.(rand(Int, length(tbd.assignment) - 5000), length(tbd.components)) .+ 1;
            tbd.assignment[1:1000] .= 0
            as1 = deepcopy(tbd.assignment);

            BPR.split_cells_by_connected_components!(tbd)
            as2 = deepcopy(tbd.assignment);

            BPR.split_cells_by_connected_components!(tbd)
            as3 = deepcopy(tbd.assignment);

            @test any(as1 .!= as2)
            @test all(sort(unique(as1)) .== sort(unique(as2)))
            @test all(sort(unique(as2)) .== sort(unique(as3)))
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
                BPR.bmm!(bm_data; n_iters=500, min_molecules_per_cell=10, assignment_history_depth=100, verbose=false);

                conj_table = counts(BPR.estimate_assignment_by_history(bm_data)[1], df_spatial.cell)
                match_masks = [vec(mapslices(maximum, conj_table ./ sum(conj_table, dims=d), dims=d)) .> 0.9 for d in 1:2]

                @test all(match_masks[2]) # We expect no new cells from multiple sources
                @test (mean(match_masks[1]) .> 0.8) # But in rare cases, we allow one cell to be split into two
                # This splitting behaviour gets rarer as we increase `n_iters`. So, we can use its frequency as a benchmark when improving the algorithm

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