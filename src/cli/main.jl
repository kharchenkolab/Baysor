using Statistics

## Configs and parameters

"""
Run cell segmentation

# Args

- `coordinates`:            CSV file with coordinates of molecules and gene type
- `prior_segmentation`:     Image or a MAT file with segmentation mask (either boolean or component indexing) or CSV
                            column with integer segmentation labels.
                            If it's the column name, it should be preceded ':' symbol (e.g. `:cell`)

# Options

- `-c, --config=<config.toml>`:         TOML file with a config
- `-x, --x-column=<x>`:                 Name of x column. Overrides the config value.
- `-y, --y-column=<y>`:                 Name of y column. Overrides the config value.
- `-z, --z-column=<z>`:                 Name of z column. Overrides the config value.
- `-g, --gene-column=<gene>`:           Name of gene column. Overrides the config value.
- `-m, --min-molecules-per-cell=<m>`:   Minimal number of molecules for a cell to be considered as real.
                                        It's an important parameter, as it's used to infer several other parameters.
                                        Overrides the config value.

- `--n-clusters=<nc>`:                  Number of molecule clusters, i.e. major cell types. Depends on protocol resolution,
                                        but should not be too high. In most cases something between 3 and 15 should work well.
                                        (default: 4)
- `--num-cells-init=<nci>`:             Initial number of cells
- `-o, --output=<path>`:                Name of the output file or path to the output directory (default: "segmentation.csv")
- `--save-polygons=<format>`:           Save estimated cell boundary polygons to a file with a specified `format`.
                                        Only 'GeoJSON' format is currently supported.
- `--scale-std=<ss>`:                   Standard deviation of scale across cells. Can be either number, which means absolute value of
                                        the std, or string ended with '%' to set it relative to scale (default: "25%")
- `-s, --scale=<s>`:                    Scale parameter, which suggest approximate cell radius for the algorithm. Must be in the same
                                        units as `x` and `y` molecule coordinates. Overrides the config value.
                                        Sets `estimate-scale-from-centers` to `false`.
- `--prior-segmentation-confidence=<p>`:    Confidence of the `prior_segmentation` results. Value in [0; 1].
                                            If you want the final segmentation not contradicting to `prior_segmentation`, set it to 1.
                                            Otherwise, if you assume errors in prior_segmentation, values in [0.2-0.7] allow
                                            flexibility for the algorithm. (default: 0.2)
- `--count-matrix-format=<format>`:     Storage format of the segmentec cell count matrix. Either 'loom' or 'tsv' (default: 'loom')


# Flags

- `-p, --plot`:             Save pdf with plot of the segmentation
- `--no-ncv-estimation`:    Turns off neighborhood composition vectors estimation

"""
@cast function run(
        coordinates::String, prior_segmentation::String="";
        config::RunOptions=RunOptions(),
        x_column::String=config.data.x, y_column::String=config.data.y, z_column::String=config.data.z,
        gene_column::String=config.data.gene, min_molecules_per_cell::Int=config.data.min_molecules_per_cell,

        scale::Float64=config.segmentation.scale, scale_std::String=config.segmentation.scale_std,
        n_clusters::Int=config.segmentation.n_clusters,
        prior_segmentation_confidence::Float64=config.segmentation.prior_segmentation_confidence,

        output::String="segmentation.csv", plot::Bool=false, save_polygons::String="false", no_ncv_estimation::Bool=false,
        count_matrix_format::String="loom"
    )

    # Parse options

    opts = config; # `config` is purely for UI clarity
    opts.data = from_dict(DataOptions,
        merge(to_dict(opts.data), Dict(
            "x" => x_column, "y" => y_column, "z" => z_column, "gene" => gene_column,
            "min_molecules_per_cell" => min_molecules_per_cell
        ))
    )

    opts.segmentation = from_dict(SegmentationOptions,
        merge(to_dict(opts.segmentation), Dict(
            "scale" => scale, "scale_std" => scale_std, "n_clusters" => n_clusters,
            "prior_segmentation_confidence" => prior_segmentation_confidence
        ))
    )

    @assert count_matrix_format in ["loom", "tsv"] "count_matrix_format must be either 'loom' or 'tsv'"

    opts, output = fill_and_check_options!(opts, prior_segmentation, output)

    # Dump run info

    Random.seed!(1)
    dump_parameters(opts, join(ARGS[2:end], " "), output)

    log_file = setup_logger(output, "log.log")

    run_id = get_run_id()
    @info "Run $run_id"
    @info get_baysor_run_str()

    # Load data

    @info "Loading data..."
    df_spatial, gene_names, n_molecules = DAT.load_df(coordinates, opts.data, filter_cols=false)

    fill_and_check_options!(opts.plotting, opts.data.min_molecules_per_cell, length(gene_names))

    if opts.segmentation.n_cells_init <= 0
        opts.segmentation.n_cells_init = default_param_value(:n_cells_init, opts.data.min_molecules_per_cell; n_molecules=n_molecules)
    end

    prior_polygons = load_prior_segmentation!(
        df_spatial, prior_segmentation, opts.segmentation; # it may change scale in opts.segmentation
        min_molecules_per_cell=opts.data.min_molecules_per_cell,
        min_molecules_per_segment=opts.data.min_molecules_per_segment, plot=plot
    )

    # Run segmentation

    @info "Estimating noise level"
    BPR.append_confidence!(
        df_spatial; nn_id=opts.data.confidence_nn_id,
        prior_confidence=opts.segmentation.prior_segmentation_confidence
    )
    @info "Done"

    (segmented_df, tracer, mol_clusts, comp_segs, poly_joint, cell_stat_df, cm, polygons) = BPR.run_segmentation(
        df_spatial, gene_names, opts.segmentation; plot_opts=opts.plotting,
        min_molecules_per_cell=opts.data.min_molecules_per_cell, estimate_ncvs=!no_ncv_estimation,
        save_polygons=(save_polygons != "false"), run_id, plot
    )

    # Save and plot results

    @info "Saving results to $output"

    out_paths = get_output_paths(output; count_matrix_format=count_matrix_format)
    DAT.save_segmentation_results(
        segmented_df, cell_stat_df, cm, polygons, out_paths;
        poly_format=save_polygons, matrix_format=Symbol(count_matrix_format), gene_names
    )

    if plot
        @info "Plotting results"
        REP.plot_segmentation_report(
            segmented_df; tracer=tracer, clust_res=mol_clusts, comp_segs=comp_segs,
            prior_polygons=prior_polygons, polygons=collect(values(poly_joint)),
            diagnostic_file=out_paths.diagnostic_report, molecule_file=out_paths.molecule_plot,
            plot_transcripts=!no_ncv_estimation, gene_colors=:ncv_color,
            min_molecules_per_cell=opts.data.min_molecules_per_cell,
            min_pixels_per_cell=opts.plotting.min_pixels_per_cell
        )
    end

    @info "All done!"

    close(log_file)

    return 0
end

function fill_and_check_options!(opts::RunOptions, prior_segmentation::String="", output::String="./")
    fill_and_check_options!(opts.data)

    if opts.segmentation.scale > 0
        opts.segmentation.estimate_scale_from_centers = false
    end

    if (prior_segmentation == "") && (opts.segmentation.scale <= 0)
        cmd_error("Either `prior_segmentation` or `scale` must be provided.")
    end

    if isdir(output) || isdirpath(output)
        output = joinpath(output, "segmentation.csv")
    end

    if xor(isempty(opts.segmentation.nuclei_genes), isempty(opts.segmentation.cyto_genes))
        cmd_error("Only one of `nuclei-genes` and `cyto-genes` is provided. It has to be either both or none.")
    end

    if (!isempty(opts.segmentation.nuclei_genes)) && (opts.segmentation.n_clusters > 1)
        @warn "Setting n-clusters > 1 is not recommended with compartment-specific expression patterns (nuclei- and cyto-genes parameters)."
    end

    return opts, output
end

function dump_parameters(options::RunOptions, args_str::String, output::String)
    open(append_suffix(output, "params.dump.toml"), "w") do f
        println(f, "# CLI params: `$args_str`")
        to_toml(f, options)
    end
end

function load_prior_segmentation!(
        df_spatial, prior_segmentation::String, opts::SegmentationOptions;
        min_molecules_per_cell::Int, min_molecules_per_segment::Int, plot::Bool
    )
    prior_polygons = Matrix{Float64}[]
    if prior_segmentation !== ""
        prior_seg_labels, scale, scale_std = DAT.load_prior_segmentation!(
            prior_segmentation, df_spatial, BPR.position_data(df_spatial);
            min_mols_per_cell=min_molecules_per_cell, min_molecules_per_segment=min_molecules_per_segment
        )

        if opts.estimate_scale_from_centers
            opts.scale, opts.scale_std = scale, string(scale_std)
        end

        if (prior_seg_labels !== nothing) && plot
            @info "Estimating prior segmentation polygons..."
            prior_polygons = BPR.boundary_polygons_from_grid(Matrix{UInt32}(prior_seg_labels[1:5:end, 1:5:end]); grid_step=5.0) # subset to save memory and time
            @info "Done"
        end
    end
    GC.gc()

    return prior_polygons
end
