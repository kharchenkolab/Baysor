import HDF5

struct OutputPaths
    cell_stats::String;
    counts::String;
    diagnostic_report::String;
    molecule_plot::String;
    polygons::String;
end
OutputPaths(;cell_stats::String, counts::String, diagnostic_report::String, molecule_plot::String, polygons::String) =
    OutputPaths(cell_stats, counts, diagnostic_report, molecule_plot, polygons)

function parse_prior_assignment(df_spatial::DataFrame, args::Dict{String, Any})
    prior_col = Symbol(args["prior_segmentation"][2:end])
    prior_segmentation = df_spatial[!, prior_col]
    try
        prior_segmentation = Int.(prior_segmentation);
    catch
        error("The prior segmentation column '$prior_col' must be of integer type")
    end
    if minimum(prior_segmentation) < 0
        error("The prior segmentation column '$prior_col' must not contain negative numbers")
    end

    filter_segmentation_labels!(prior_segmentation, min_molecules_per_segment=args["min-molecules-per-segment"])
    if args["scale"] === nothing
        pos_data = position_data(df_spatial)
        scale, scale_std = estimate_scale_from_assignment(pos_data, prior_segmentation;
            min_mols_per_cell=args["min-molecules-per-cell"])
    else
        scale, scale_std = args["scale"], args["scale-std"]
    end

    return prior_segmentation, Matrix{Float64}[], scale, scale_std
end

function load_prior_segmentation(df_spatial::DataFrame, args::Dict{String, Any})
    if args["prior_segmentation"][1] == ':'
        return parse_prior_assignment(df_spatial, args)
    end

    @info "Loading segmentation mask..."

    prior_seg_labels = load_segmentation_mask(args["prior_segmentation"])
    prior_segmentation = Int.(staining_value_per_transcript(df_spatial, prior_seg_labels));
    filter_segmentation_labels!(prior_seg_labels, prior_segmentation; min_molecules_per_segment=args["min-molecules-per-segment"])
    GC.gc()

    scale, scale_std = args["scale"], args["scale-std"]
    if args["estimate-scale-from-centers"]
        scale, scale_std = estimate_scale_from_centers(prior_seg_labels)
    end
    @info "Done"

    prior_polygons = nothing
    if args["plot"]
        @info "Estimating prior segmentation polygons..."
        prior_polygons = extract_polygons_from_label_grid(Matrix{UInt32}(prior_seg_labels[1:5:end, 1:5:end]); grid_step=5.0) # subset to save memory and time
        @info "Done"
    end

    return prior_segmentation, prior_polygons, scale, scale_std
end

### Data

function save_segmentation_results(bm_data::BmmData, gene_names::Vector{String}, args::Dict{String, Any}; paths::OutputPaths,
        mol_clusts::Union{NamedTuple, Nothing}, comp_segs::Union{NamedTuple, Nothing}, prior_polygons::Union{Vector{Matrix{Float64}}, Nothing})
    segmentated_df = get_segmentation_df(bm_data, gene_names)
    cell_stat_df = get_cell_stat_df(bm_data, segmentated_df; add_qc=true)

    if args["estimate-ncvs"]
        @info "Estimating local colors"
        gene_colors = gene_composition_colors(bm_data.x, args["gene-composition-neigborhood"])
        segmentated_df[!, :ncv_color] = "#" .* Colors.hex.(gene_colors)
    end

    @info "Saving results to $(args["output"])"
    CSV.write(args["output"], segmentated_df[sortperm(segmentated_df.molecule_id), :]);
    CSV.write(paths.cell_stats, cell_stat_df);

    cm = convert_segmentation_to_counts(bm_data.x.gene, bm_data.assignment; gene_names=gene_names)
    count_str = join(names(cm), "\t") * "\n" * join([join(["$v" for v in r], '\t') for r in eachrow(cm)], '\n');
    open(paths.counts, "w") do f; print(f, count_str) end

    if args["plot"]
        REP.plot_diagnostics_panel(
            segmentated_df, segmentated_df.cell, bm_data.tracer;
            file=paths.diagnostic_report,
            clust_res=mol_clusts, comp_segs=comp_segs
        )
        if args["estimate-ncvs"]
            poly_joint, polygons = boundary_polygons_auto(
                bm_data.x, bm_data.assignment; scale=args["scale"], min_pixels_per_cell=args["min-pixels-per-cell"],
                estimate_per_z=(args["save-polygons"] !== nothing)
            )
            REP.plot_transcript_assignment_panel(
                bm_data.x; clusters=bm_data.cluster_per_molecule, gene_colors=gene_colors, file=paths.molecule_plot,
                min_molecules_per_cell=args["min-molecules-per-cell"], min_pixels_per_cell=args["min-pixels-per-cell"],
                prior_polygons=prior_polygons, polygons=poly_joint,
            )

            if args["save-polygons"] !== nothing
                if lowercase(args["save-polygons"]) == "geojson"
                    save_polygons_to_geojson(polygons, paths.polygons)
                else
                    @warn "Unknown polygon format: $(args["save-polygons"])"
                end
            end
        end
    end
end

function save_matrix_to_loom(matrix; gene_names::Vector{String}, cell_names::Vector{String}, file_path::String,
        row_attrs::Union{Dict{String, T1}, Nothing} where T1=nothing, col_attrs::Union{Dict{String, T2}, Nothing} where T2=nothing)
    # Specification: https://linnarssonlab.org/loompy/format/index.html
    HDF5.h5open(file_path, "w") do fid
        fid["matrix", chunk=(64,64), compress=3] = matrix
        HDF5.create_group(fid, "row_attrs")
        HDF5.create_group(fid, "col_attrs")
        HDF5.create_group(fid, "attrs")
        fid["row_attrs"]["Name"] = gene_names
        fid["col_attrs"]["CellID"] = cell_names
        if row_attrs !== nothing
            for (k,v) in row_attrs
                fid["row_attrs"][k] = v
            end
        end

        if col_attrs !== nothing
            for (k,v) in col_attrs
                fid["col_attrs"][k] = v
            end
        end
    end;
end

### Processing

function estimate_molecule_clusters(df_spatial::DataFrame, n_clusters::Int)
    @info "Clustering molecules..."
    # , adjacency_type=:both, k_adj=fmax(1, div(args["min-molecules-per-cell"], 2))
    adjacent_points, adjacent_weights = build_molecule_graph_normalized(df_spatial, :confidence, filter=false);

    mol_clusts = cluster_molecules_on_mrf(df_spatial, adjacent_points, adjacent_weights; n_clusters=n_clusters, weights_pre_adjusted=false)

    @info "Done"
    return mol_clusts
end

function estimate_molecule_compartments(df_spatial::DataFrame, gene_names::Vector{String}, args::Dict{String, Any})
    @info "Estimating compartment regions..."
    comp_genes = Dict{String, Vector{String}}()

    # Parse genes from CLI
    for k in ["nuclei-genes", "cyto-genes"]
        c_genes = String.(strip.(Base.split(args[k], ",")))
        all(length.(c_genes) .> 0) || @warn "Empty genes were provided in $k"

        missing_genes = c_genes[.!in.(c_genes, Ref(Set(gene_names)))]
        (length(missing_genes) == 0) || @warn "Genes $(join(missing_genes, ',')) are missing from the data"
        c_genes = intersect(c_genes, gene_names)
        length(c_genes) > 0 || error("No genes left in $k after filtration")
        comp_genes[k] = c_genes
    end

    # Run segmentation
    adjacent_points, adjacent_weights = build_molecule_graph_normalized(df_spatial, :confidence, filter=false);

    init_probs, is_locked = init_nuclei_cyto_compartments(
        position_data(df_spatial), df_spatial.gene; gene_names=gene_names, scale=args["scale"],
        nuclei_genes=comp_genes["nuclei-genes"], cyto_genes=comp_genes["cyto-genes"]
    );

    comp_segs = segment_molecule_compartments(init_probs, is_locked, adjacent_points, adjacent_weights, df_spatial.confidence);

    @info "Done"

    id_per_gene = Dict(g => i for (i,g) in enumerate(gene_names))
    comp_genes = [[id_per_gene[g] for g in gs] for gs in values(comp_genes)]

    return comp_segs, comp_genes
end
