function append_confidence!(df_spatial::DataFrame, args::Dict{String})
    confidence_nn_id = default_param_value(:confidence_nn_id, args["min-molecules-per-cell"])

    @info "Estimating noise level"
    prior_seg = (args["prior_segmentation"]===nothing) ? nothing : df_spatial.prior_segmentation
    append_confidence!(df_spatial, prior_seg, nn_id=confidence_nn_id, prior_confidence=args["prior-segmentation-confidence"])
    @info "Done"
end

function estimate_gene_structure_embedding(df_spatial::DataFrame, gene_names::Vector{String}, confidence::Vector{Float64}=df_spatial.confidence)
    adjacent_points, adjacent_weights = build_molecule_graph(df_spatial, filter=false)[1:2];
    cor_mat = pairwise_gene_spatial_cor(df_spatial.gene, confidence, adjacent_points, adjacent_weights);
    cor_vals = vec(cor_mat)
    min_cor = quantile(cor_vals[cor_vals .> 0], 0.001) ./ 5;
    max_cor = quantile(cor_vals, 0.99);
    p_dists = 1 .- max.(min.(cor_mat, max_cor), min_cor) ./ max_cor;
    p_dists[diagind(p_dists)] .= 0.0;

    embedding = UMAP.umap(p_dists, 2; metric=:precomputed, spread=1.0, min_dist=0.1, n_epochs=5000, n_neighbors=max(min(15, length(gene_names) ÷ 2), 2));
    marker_sizes = log.(count_array(df_spatial.gene));

    return DataFrame(Dict(:x => embedding[1,:], :y => embedding[2,:], :gene => Symbol.(gene_names), :size => marker_sizes));
end

function run_segmentation(df_spatial::DataFrame, gene_names::Vector{String}, args::Dict{String, Any})
    # run_id = get_run_id()
    # @info "Run $run_id"
    # TODO: add run_id to cell ids

    # Region-based segmentations

    comp_segs, comp_genes = nothing, Vector{Int}[]
    adjacent_points, adjacent_weights = build_molecule_graph(df_spatial; use_local_gene_similarities=false, adjacency_type=:triangulation)[1:2]
    if args["nuclei-genes"] !== nothing
        comp_segs, comp_genes, df_spatial[!, :compartment] = estimate_molecule_compartments(
            df_spatial, gene_names; nuclei_genes=args["nuclei-genes"], cyto_genes=args["cyto-genes"], scale=args["scale"]
        )
        df_spatial[!, :nuclei_probs] = 1 .- comp_segs.assignment_probs[2,:];

        adjacent_points, adjacent_weights = adjust_mrf_with_compartments(
            adjacent_points, adjacent_weights,
            comp_segs.assignment_probs[1,:], comp_segs.assignment_probs[2,:]
        );
    end

    mol_clusts = nothing
    if args["n-clusters"] > 1
        mol_clusts = estimate_molecule_clusters(df_spatial, args["n-clusters"])
        df_spatial[!, :cluster] = mol_clusts.assignment;
    end

    # Cell segmentation

    n_cells_init = something(
        args["num-cells-init"],
        default_param_value(:n_cells_init, args["min-molecules-per-cell"], n_molecules=size(df_spatial, 1))
    )

    bm_data = initialize_bmm_data(
        df_spatial; scale=args["scale"], scale_std=args["scale-std"], n_cells_init=n_cells_init,
        prior_seg_confidence=args["prior-segmentation-confidence"], min_molecules_per_cell=args["min-molecules-per-cell"],
        adjacent_points=adjacent_points, adjacent_weights=adjacent_weights, na_genes=Vector{Int}(vcat(comp_genes...))
    );

    @info "Using $(size(position_data(bm_data), 1))D coordinates"

    history_depth = round(Int, args["iters"] * 0.1)
    bm_data = bmm!(
        bm_data; n_iters=args["iters"],
        new_component_frac=args["new-component-fraction"], new_component_weight=args["new-component-weight"],
        min_molecules_per_cell=args["min-molecules-per-cell"], assignment_history_depth=history_depth
    );

    @info "Processing complete."

    # Extract results

    segmented_df, cell_stat_df, cm = get_segmentation_results(bm_data, gene_names)
    gene_colors = nothing
    if args["estimate-ncvs"]
        @info "Estimating local colors"
        gene_colors = gene_composition_colors(bm_data.x, args["gene-composition-neigborhood"])
        segmented_df[!, :ncv_color] = "#" .* Colors.hex.(gene_colors)
    end

    poly_joint, polygons = nothing, nothing
    if args["save-polygons"] !== nothing || args["plot"] && args["estimate-ncvs"]
        poly_joint, polygons = boundary_polygons_auto(
            position_data(bm_data), bm_data.assignment; scale=args["scale"], min_pixels_per_cell=args["min-pixels-per-cell"],
            estimate_per_z=(args["save-polygons"] !== nothing)
        )
    end

    return (
        segmented_df, bm_data.tracer, mol_clusts, comp_segs, poly_joint,
        cell_stat_df, cm, polygons
    )
end