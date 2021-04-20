using ArgParse
using DataFrames
using DataFramesMeta
using ProgressMeter
using Statistics

import Colors
import CSV
import Dates
import LibGit2
import Pkg
import Pkg.TOML
import Plots

## Common

parse_toml_config(config::T where T <: AbstractString) =
    parse_toml_config(TOML.parsefile(config))

function get_default_config()
    return deepcopy(Dict{String, Any}(
        "Data" => Dict{String, Any}(
            "x-column" => "x",
            "y-column" => "y",
            "z-column" => "z",
            "gene-column" => "gene",
            "min-molecules-per-gene" => 1,
            "min-molecules-per-cell" => 3,
            "estimate-scale-from-centers" => true,
            "scale" => nothing,
            "scale-std" => "25%",
            "min-molecules-per-segment" => nothing
        ),
        "Sampling" => Dict{String, Any}(
            "new-component-weight" => 0.2,
            "new-component-fraction" => 0.3
        ),
        "Plotting" => Dict{String, Any}(
            "gene-composition-neigborhood" => nothing,
            "min-pixels-per-cell" => 15
        )
    ))
end

function parse_toml_config(config::Dict{AS, Any}) where AS <: AbstractString
    res_config = get_default_config()
    for (k,v) in config
        if !(k in keys(res_config))
            error("Unexpected value in the config: '$k'")
        end

        cur_def = res_config[k]

        for (k2,v2) in v
            if !(k2 in keys(cur_def))
                error("Unexpected value in the config: '$k' -> '$k2'")
            end

            cur_def[k2] = v2
        end
    end

    return res_config
end

function extend_params_with_config!(params::Dict, config::Dict)
    for sub_cfg in values(config)
        for (k, v) in sub_cfg
            if !(k in keys(params)) || params[k] === nothing
                params[k] = v
            end
        end
    end
end

function load_df(args::Dict; kwargs...) 
    exc_genes = (args["exclude-genes"] === nothing) ? String[] : String.(strip.(Base.split(args["exclude-genes"], ",")))

    df_spatial, gene_names = load_df(args["coordinates"]; x_col=args["x-column"], y_col=args["y-column"], z_col=args["z-column"], gene_col=args["gene-column"], 
        min_molecules_per_gene=args["min-molecules-per-gene"], exclude_genes=exc_genes, kwargs...)
    
    if args["force-2d"] && (:z in propertynames(df_spatial))
        select!(df_spatial, Not(:z))
    elseif (args["z-column"] != :z) && !(:z in propertynames(df_spatial))
        error("z-column $(args["z-column"]) not found in the data")
    end

    return df_spatial, gene_names
end

append_suffix(output::String, suffix) = "$(splitext(output)[1])_$suffix"

## Main run

function parse_commandline(args::Union{Nothing, Array{String, 1}}=nothing) # TODO: add verbosity level
    s = ArgParseSettings(prog="baysor run")
    @add_arg_table! s begin
        "--config", "-c"
            help = "TOML file with config"
        "--x-column", "-x"
            help = "Name of x column. Overrides the config value."
        "--y-column", "-y"
            help = "Name of y column. Overrides the config value."
        "--z-column", "-z"
            help = "Name of z column. Overrides the config value."
        "--gene-column", "-g"
            help = "Name of gene column. Overrides the config value."

        "--iters", "-i"
            help = "Number of iterations"
            arg_type = Int
            default = 500
        "--min-molecules-per-cell", "-m"
            help = "Minimal number of molecules for a cell to be considered as real. It's an important parameter, as it's used to infer several other parameters. Overrides the config value."
            arg_type = Int
        "--n-clusters"
            help = "Number of molecule clusters, i.e. major cell types. Depends on protocol resolution, but should not be too high. In most cases something between 3 and 15 should work well."
            arg_type = Int
            default=4
        "--num-cells-init"
            help = "Initial number of cells."
            arg_type = Int
        "--output", "-o"
            help = "Name of the output file or path to the output directory"
            default = "segmentation.csv"
        "--plot", "-p"
            help = "Save pdf with plot of the segmentation"
            action = :store_true
        "--save-polygons"
            help = "Save estimated cell boundary polygons to a file with a specified FORMAT. Only 'GeoJSON' format is currently supported. The option requires setting '-p' to work."
            range_tester = (x -> in(lowercase(x), ["geojson"]))
            metavar = "FORMAT"
        "--force-2d"
            help = "Ignores z-column in the data if it is provided"
            action = :store_true
        "--exclude-genes"
            help = "Comma-separated list of genes to ignore during segmentation"
        "--nuclei-genes"
            help = "Comma-separated list of nuclei-specific genes. If provided, `cyto-genes` has to be set, as well."
        "--cyto-genes"
            help = "Comma-separated list of cytoplasm-specific genes. If provided, `nuclei-genes` has to be set, as well."
        "--scale-std"
            help = "Standard deviation of scale across cells. Can be either number, which means absolute value of the std, or string ended with '%' to set it relative to scale. Default: 25%"
            default = "25%"

        "--scale", "-s"
            help = "Scale parameter, which suggest approximate cell radius for the algorithm. Must be in the same units as 'x' and 'y' molecule coordinates. Overrides the config value. Sets 'estimate-scale-from-centers' to false."
            arg_type = Float64
        "--prior-segmentation-confidence"
            help = "Confidence of the `prior_segmentation` results. Value in [0; 1]. If you want the final segmentation not contradicting to prior_segmentation, set it to 1. Otherwise, if you assume errors in prior_segmentation, values in [0.2-0.7] allow flexibility for the algorithm."
            arg_type = Float64
            default = 0.2

        "coordinates"
            help = "CSV file with coordinates of molecules and gene type"
            required = true
        "prior_segmentation"
            help = "Image or a MAT file with segmentation mask (either boolean or component indexing) or CSV column with integer segmentation labels. If it's the column name, it should be preceded ':' symbol (e.g. :cell)"
    end

    return (args === nothing) ? parse_args(s) : parse_args(args, s)
end

function parse_configs(args::Union{Nothing, Array{String, 1}}=nothing)
    r = parse_commandline(args)
    if r["scale"] !== nothing
        r["estimate-scale-from-centers"] = false
    end

    if r["config"] !== nothing
        extend_params_with_config!(r, parse_toml_config(r["config"]))
    else
        @warn "No config file provided. Back-up to default parameters."
        extend_params_with_config!(r, get_default_config())
    end

    for k in ["gene-column", "x-column", "y-column"]
        if !(k in keys(r)) || (r[k] === nothing) # should never be the case as we have defaults
            @warn "$k must be specified"
            return nothing
        end
        r[k] = Symbol(r[k])
    end

    r["z-column"] = Symbol(r["z-column"])

    if r["prior_segmentation"] === nothing && r["scale"] === nothing
        @warn "Either `prior_segmentation` or `scale` must be provided."
        return nothing
    end

    if r["min-molecules-per-segment"] === nothing
        r["min-molecules-per-segment"] = default_param_value(:min_molecules_per_segment, r["min-molecules-per-cell"])
    end

    if isdir(r["output"]) || isdirpath(r["output"])
        r["output"] = joinpath(r["output"], "segmentation.csv")
    end

    if (r["save-polygons"] !== nothing) && (!r["plot"])
        @warn "--plot option is required for saving polygons (--save-polygons). The polygons will not be saved."
    end

    if xor(r["nuclei-genes"] === nothing, r["cyto-genes"] === nothing)
        @warn "Only one of `nuclei-genes` and `cyto-genes` is provided. It has to be either both or none."
        return nothing
    end

    if (r["nuclei-genes"] !== nothing) && (r["n-clusters"] > 1)
        @warn "Setting n-clusters > 1 is not recommended with compartment-specific expression patterns (nuclei- and cyto-genes parameters)."
    end

    return r
end

function plot_diagnostics_panel(df_res::DataFrame, assignment::Array{Int, 1}, tracer::Dict, args::Dict; margin=5*Plots.mm, 
        clust_res::Union{NamedTuple, Nothing}=nothing, comp_segs::Union{NamedTuple, Nothing}=nothing)
    @info "Plot diagnostics"
    open(append_suffix(args["output"], "diagnostics.html"), "w") do io
        vega_plots = Dict{String, VegaLite.VLSpec}()
        println(io, "<html>")
        println(io, vega_header("Report"))
        println(io, vega_style())
        println(io, "<body>")
        # Molecule clustering convergence
        if clust_res !== nothing
            println(io, "<div id='vg_clust_conv'></div>")
            vega_plots["vg_clust_conv"] = plot_clustering_convergence(clust_res, "Molecule clustering convergence")
        end

        if comp_segs !== nothing
            println(io, "<div id='vg_compart_conv'></div>")
            vega_plots["vg_compart_conv"] = plot_clustering_convergence(comp_segs, "Compartment segmentation convergence")
        end

        println(io, "<br><br>")
        # Main algorithm convergence
        if (:n_components in keys(tracer)) && length(tracer[:n_components]) != 0
            p_conv = plot_num_of_cells_per_iterarion(tracer);
            show(io, MIME("text/html"), p_conv)
        end

        println(io, "<br><br>")

        # Confidence per molecule
        if :confidence in propertynames(df_res)
            bins = 0.0:0.025:1.0
            p_conf = Plots.histogram(df_res.confidence[assignment .!= 0], bins=bins, label="Assigned molecules",
                xlabel="Confidence", ylabel="#Molecules", title="Confidence per molecule", margin=margin)
            p_conf = Plots.histogram!(df_res.confidence[assignment .== 0], alpha=0.5, bins=bins, label="Noise molecules")
            show(io, MIME("text/html"), p_conf)
        end

        # Assignment confidence
        if :assignment_confidence in propertynames(df_res)
            p_conf = Plots.histogram(df_res.assignment_confidence[assignment .> 0], bins=50, legend=false)
            p_conf = Plots.vline!([0.95], xlabel="Assignment confidence", ylabel="#Molecules", xlims=(-0.01, 1.03),
                title="Assignment confidence per real molecules")
            show(io, MIME("text/html"), p_conf)
        end

        println(io, "<br><br>")

        # Num. of molecules per cell
        n_mols_per_cell = count_array(assignment, drop_zero=true)
        p_n_mols = Plots.histogram(n_mols_per_cell[(n_mols_per_cell .> 1) .& (n_mols_per_cell .< quantile(n_mols_per_cell, 0.99) / 0.99)],
            title="Num. molecules per cell", xlabel="Num. molecules per cell", ylabel="Num. cells", label=:none)
        show(io, MIME("text/html"), p_n_mols)

        println(io, "</body>")
        if length(vega_plots) > 0
            println(io, vega_style())
            println(io, vega_plot_html(vega_plots))
        end
        println(io, "</html>")
    end
end

function plot_transcript_assignment_panel(df_res::DataFrame, assignment::Vector{Int}, args::Dict;
        clusters::Vector{Int}, prior_polygons::Array{Matrix{Float64}, 1}, gene_colors::Union{Vector{<:Colors.ColorTypes.Color}, Nothing}=nothing)
    if gene_colors === nothing
        @info "Estimating local colors"
        gene_colors = gene_composition_colors(df_res, args["gene-composition-neigborhood"])
    end

    @info "Estimating boundary polygons" # For some parameters, this step takes a lot of time and memory
    grid_step = args["scale"] / args["min-pixels-per-cell"];
    polygons = boundary_polygons(df_res, assignment; grid_step=grid_step, bandwidth=args["scale"]/10);

    if gene_colors !== nothing
        gene_colors = Colors.alphacolor.(gene_colors, 0.5)
    end

    @info "Plot transcript assignment"
    gc_plot = plot_dataset_colors(df_res, gene_colors; polygons=polygons, prior_polygons=prior_polygons, min_molecules_per_cell=args["min-molecules-per-cell"],
        min_pixels_per_cell=args["min-pixels-per-cell"], title="Local expression similarity")

    clust_plot = nothing
    if !isempty(clusters)
        clust_plot = plot_dataset_colors(df_res, gene_colors; polygons=polygons, prior_polygons=prior_polygons, annotation=clusters, min_molecules_per_cell=args["min-molecules-per-cell"],
            min_pixels_per_cell=args["min-pixels-per-cell"], title="Molecule clustering")
    end

    open(append_suffix(args["output"], "borders.html"), "w") do io
        show(io, MIME("text/html"), gc_plot)

        if clust_plot !== nothing
            show(io, MIME("text/html"), clust_plot)
        end
    end

    return polygons
end

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

    @info "Estimating prior segmentation polygons..."
    prior_polygons = extract_polygons_from_label_grid(Matrix{UInt32}(prior_seg_labels[1:3:end, 1:3:end]); grid_step=3.0) # subset to save memory and time
    @info "Done"

    return prior_segmentation, prior_polygons, scale, scale_std
end

function get_baysor_run_str()::String
    pkg_info = Pkg.dependencies()[Base.UUID("cc9f9468-1fbe-11e9-0acf-e9460511877c")]
    pkg_str = "v$(pkg_info.version)"
    try
        repo = LibGit2.GitRepo(pkg_info.source)
        hash_str = "$(LibGit2.GitShortHash(LibGit2.peel(LibGit2.GitCommit, LibGit2.head(repo))))"
        pkg_str = "$pkg_str [$hash_str]"
    catch err
    end

    return "($(Dates.Date(Dates.now()))) Run Baysor $pkg_str"
end

function dump_parameters(args::Dict{String, Any}, args_str::String)
    if args["config"] !== nothing
        dump_dst = abspath(append_suffix(args["output"], "config.toml"))
        if abspath(args["config"]) != dump_dst
            cp(args["config"], dump_dst, force=true)
        end
    end

    open(append_suffix(args["output"], "params.dump"), "w") do f
        println(f, "# CLI params: `$args_str`")
        TOML.print(f, Dict(k => (v !== nothing) ? ((typeof(v) === Symbol) ? String(v) : v) : "" for (k,v) in args))
    end
end

function load_and_preprocess_data!(args::Dict{String, Any})
    @info get_baysor_run_str()
    @info "Loading data..."
    df_spatial, gene_names = load_df(args, filter_cols=false)
    df_spatial[!, :molecule_id] = 1:size(df_spatial, 1)

    @info "Loaded $(size(df_spatial, 1)) transcripts"
    if size(df_spatial, 1) != size(unique(df_spatial), 1)
        @warn "$(size(df_spatial, 1) - size(unique(df_spatial), 1)) records are duplicates. You may need to filter them beforehand."
    end

    if args["gene-composition-neigborhood"] === nothing
        args["gene-composition-neigborhood"] = default_param_value(:composition_neighborhood, args["min-molecules-per-cell"], n_genes=length(gene_names))
    end

    prior_polygons = Matrix{Float64}[]
    if args["prior_segmentation"] !== nothing
        df_spatial[!, :prior_segmentation], prior_polygons, args["scale"], args["scale-std"] = load_prior_segmentation(df_spatial, args)
    end
    GC.gc()

    confidence_nn_id = default_param_value(:confidence_nn_id, args["min-molecules-per-cell"])

    @info "Estimating noise level"
    prior_seg = (args["prior_segmentation"]===nothing) ? nothing : df_spatial.prior_segmentation
    append_confidence!(df_spatial, prior_seg, nn_id=confidence_nn_id, prior_confidence=args["prior-segmentation-confidence"])
    @info "Done"

    return df_spatial, gene_names, prior_polygons
end

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
    comp_genes = Vector{String}[]

    # Parse genes from CLI
    for k in ["nuclei-genes", "cyto-genes"]
        c_genes = String.(strip.(Base.split(args[k], ",")))
        all(length.(c_genes) .> 0) || @warn "Empty genes were provided in $k"
        
        missing_genes = c_genes[.!in.(c_genes, Ref(Set(gene_names)))]
        (length(missing_genes) == 0) || @warn "Genes $(join(missing_genes, ',')) are missing from the data"
        c_genes = intersect(c_genes, gene_names)
        length(c_genes) > 0 || error("No genes left in $k after filtration")
        push!(comp_genes, c_genes)
    end

    # Run segmentation
    nn_id = default_param_value(:compartment_nn_id, args["min-molecules-per-cell"])
    adjacent_points, adjacent_weights = build_molecule_graph_normalized(df_spatial, :confidence, filter=false);
    comp_segs = segment_molecule_compartments(position_data(df_spatial), df_spatial.gene, 
        adjacent_points, adjacent_weights, df_spatial.confidence; comp_genes=comp_genes, gene_names=gene_names, nn_id=nn_id);

    @info "Done"

    id_per_gene = Dict(g => i for (i,g) in enumerate(gene_names))
    comp_genes = [[id_per_gene[g] for g in gs] for gs in comp_genes]

    return comp_segs, comp_genes
end

function save_segmentation_results(bm_data::BmmData, gene_names::Vector{String}, args::Dict{String, Any}; 
        mol_clusts::Union{NamedTuple, Nothing}, comp_segs::Union{NamedTuple, Nothing}, prior_polygons::Vector{Matrix{Float64}})
    segmentated_df = get_segmentation_df(bm_data, gene_names)
    cell_stat_df = get_cell_stat_df(bm_data, segmentated_df; add_qc=true)

    @info "Estimating local colors"
    gene_colors = gene_composition_colors(bm_data.x, args["gene-composition-neigborhood"])
    segmentated_df[!, :ncv_color] = "#" .* Colors.hex.(gene_colors)

    @info "Saving results to $(args["output"])"
    CSV.write(args["output"], segmentated_df[sortperm(segmentated_df.molecule_id), :]);
    CSV.write(append_suffix(args["output"], "cell_stats.csv"), cell_stat_df);

    cm = convert_segmentation_to_counts(bm_data.x.gene, bm_data.assignment; gene_names=gene_names)
    count_str = join(names(cm), "\t") * "\n" * join([join(["$v" for v in r], '\t') for r in eachrow(cm)], '\n');
    open(append_suffix(args["output"], "counts.tsv"), "w") do f; print(f, count_str) end

    if args["plot"]
        plot_diagnostics_panel(segmentated_df, segmentated_df.cell, bm_data.tracer, args; clust_res=mol_clusts, comp_segs=comp_segs)
        polygons = plot_transcript_assignment_panel(bm_data.x, bm_data.assignment, args; clusters=bm_data.cluster_per_molecule, prior_polygons=prior_polygons,
            gene_colors=gene_colors)

        if args["save-polygons"] !== nothing
            if lowercase(args["save-polygons"]) == "geojson"
                save_polygons_to_geojson(polygons, append_suffix(args["output"], "polygons.json"))
            else
                @warn "Unknown polygon format: $(args["save-polygons"])"
            end
        end
    end
end

function run_cli_main(args::Union{Nothing, Array{String, 1}}=nothing)
    args_str = join(args, " ")
    args = parse_configs(args)
    (args !== nothing) || return 1

    dump_parameters(args, args_str)

    # Set up logger

    log_file = open(append_suffix(args["output"], "log.log"), "w")
    Base.CoreLogging.global_logger(DoubleLogger(log_file, stdout; force_flush=true))

    # Set up plotting

    ENV["GKSwstype"] = "100"; # Disable output device
    Plots.gr()

    df_spatial, gene_names, prior_polygons = load_and_preprocess_data!(args)

    # Region-based segmentations

    comp_segs, comp_genes = nothing, Vector{Int}[]
    adjacent_points, adjacent_weights = build_molecule_graph(df_spatial; use_local_gene_similarities=false, adjacency_type=:triangulation)[1:2]
    if args["nuclei-genes"] !== nothing
        comp_segs, comp_genes = estimate_molecule_compartments(df_spatial, gene_names, args)
        df_spatial[!, :compartment] = ["Nuclei", "Cyto", "Unknown"][comp_segs.assignment];
        df_spatial[!, :nuclei_probs] = 1 .- comp_segs.assignment_probs[2,:];

        adjacent_points, adjacent_weights = adjust_mrf_with_compartments(adjacent_points, adjacent_weights, 
            comp_segs.assignment_probs[1,:], comp_segs.assignment_probs[2,:]);
    end

    mol_clusts = nothing
    if args["n-clusters"] > 1
        mol_clusts = estimate_molecule_clusters(df_spatial, args["n-clusters"])
        df_spatial[!, :cluster] = mol_clusts.assignment;
    end

    # Cell segmentation

    bm_data = initialize_bmm_data(df_spatial; scale=args["scale"], scale_std=args["scale-std"],
            n_cells_init=args["num-cells-init"], prior_seg_confidence=args["prior-segmentation-confidence"],
            min_molecules_per_cell=args["min-molecules-per-cell"], confidence_nn_id=0,
            adjacent_points=adjacent_points, adjacent_weights=adjacent_weights, na_genes=vcat(comp_genes...));
        
    @info "Using $(size(position_data(bm_data), 1))D coordinates"

    history_depth = round(Int, args["iters"] * 0.1)
    bm_data = bmm!(bm_data; n_iters=args["iters"], new_component_frac=args["new-component-fraction"], new_component_weight=args["new-component-weight"],
            min_molecules_per_cell=args["min-molecules-per-cell"], assignment_history_depth=history_depth);

    @info "Processing complete."

    # Save results
    save_segmentation_results(bm_data, gene_names, args; mol_clusts=mol_clusts, comp_segs=comp_segs, prior_polygons=prior_polygons)
    @info "All done!"

    close(log_file)

    return 0
end

## Pre-view

function parse_preview_commandline(args::Union{Nothing, Array{String, 1}}=nothing)
    s = ArgParseSettings(prog="baysor preview")
    @add_arg_table! s begin
        "--config", "-c"
            help = "TOML file with config"
        "--x-column", "-x"
            help = "Name of x column. Overrides the config value."
        "--y-column", "-y"
            help = "Name of gene column. Overrides the config value."
        "--gene-column", "-g"
            help = "Name of gene column. Overrides the config value."
        "--min-molecules-per-cell", "-m"
            help = "Minimal number of molecules for a cell to be considered as real. It's an important parameter, as it's used to infer several other parameters. Overrides the config value."
            arg_type = Int
        "--min-pixels-per-cell"
            help = "Minimal number of pixels per cell. Used to estimate size of the dataset plot."
            arg_type = Int
            default = 15
        "--output", "-o"
            help = "Name of the output file or path to the output directory"
            default = "preview.html"
        "coordinates"
            help = "CSV file with coordinates of transcripts and gene type"
            required = true
        # TODO: add gene-composition-neighborhood
        # TODO: add verbosity level
    end

    return (args === nothing) ? parse_args(s) : parse_args(args, s)
end

function parse_preview_configs(args::Union{Nothing, Array{String, 1}}=nothing)
    r = parse_preview_commandline(args)
    if r["config"] !== nothing
        extend_params_with_config!(r, parse_toml_config(r["config"]))
    else
        extend_params_with_config!(r, get_default_config())
    end

    for k in ["gene-column", "x-column", "y-column"]
        r[k] = Symbol(r[k])
    end

    if isdir(r["output"]) || isdirpath(r["output"])
        r["output"] = joinpath(r["output"], "preview.html")
    end

    return r
end

function run_cli_preview(args::Union{Nothing, Array{String, 1}}=nothing)
    args = parse_preview_configs(args)

    # Set up logger

    log_file = open(append_suffix(args["output"], "preview_log.log"), "w")
    Base.CoreLogging.global_logger(DoubleLogger(log_file, stdout; force_flush=true))

    # Set up plotting

    ENV["GKSwstype"] = "100"; # Disable output device
    Plots.gr()

    # Run preview

    @info get_baysor_run_str()
    @info "Loading data..."
    args["min-molecules-per-gene"] = 0
    df_spatial, gene_names = load_df(args)

    @info "Loaded $(size(df_spatial, 1)) transcripts"
    if size(df_spatial, 1) != size(unique(df_spatial), 1)
        @warn "$(size(df_spatial, 1) - size(unique(df_spatial), 1)) records are duplicates. You may need to filter them beforehand."
    end

    @info "Estimating noise level"
    confidence_nn_id = default_param_value(:confidence_nn_id, args["min-molecules-per-cell"])
    edge_lengths, confidences, d1, d2 = append_confidence!(df_spatial, nn_id=confidence_nn_id) # TODO: use segmentation mask if available here
    @info "Done"

    @info "Estimating local neighborhoods"

    if args["gene-composition-neigborhood"] === nothing
        args["gene-composition-neigborhood"] = default_param_value(:composition_neighborhood, args["min-molecules-per-cell"], n_genes=length(gene_names))
    end

    @info "Estimating local colors"
    gene_colors = gene_composition_colors(df_spatial, args["gene-composition-neigborhood"])

    @info "Building transcript plots"
    gc_plot = plot_dataset_colors(df_spatial, gene_colors; min_molecules_per_cell=args["min-molecules-per-cell"],
        min_pixels_per_cell=args["min-pixels-per-cell"], title="Local expression similarity")

    conf_colors = map_to_colors(confidences, lims=(0.0, 1.0), palette=Colors.diverging_palette(10, 250, s=0.75, w=1.0));
    cc_plot = plot_dataset_colors(df_spatial, conf_colors[:colors]; min_molecules_per_cell=args["min-molecules-per-cell"],
        min_pixels_per_cell=args["min-pixels-per-cell"], title="Transcript confidence")

    @info "Building gene structure plot"
    gene_structure_plot = plot_gene_structure(df_spatial, gene_names, format=:png)
    ## Plots

    n_tr_plot = plot_num_transcript_overview(df_spatial.gene, confidences, gene_names)
    noise_dist_plot = plot_noise_estimation_diagnostics(edge_lengths, confidences, d1, d2, confidence_nn_id=confidence_nn_id)

    @info "Plotting"

    open(args["output"], "w") do io
        print(io, """
            <!DOCTYPE html>
            <html>
            <head><style> body {background-color: #ffffff;}</style></head>

            <body>
            <div>
            <h2>Content</h2>
            <ul>
                <li><a href="#Transcript_Plots">Transcript plots</a></li>
                <li><a href="#Noise_Level">Noise level</a></li>
                <li><a href="#Gene_Structure">Gene structure</a></li>
            </ul>
            </div>
            """)

        println(io, "<h1 id=\"Transcript_Plots\">Transcript plots</h1><br>")
        show(io, MIME("text/html"), gc_plot)
        println(io, "<br>")

        println(io, "<hr>\n<h1 id=\"Noise_Level\">Noise level</h1><br>")
        show(io, MIME("text/html"), cc_plot)
        println(io, "<br>")
        show(io, MIME("text/html"), noise_dist_plot)
        println(io, "<br>")
        println(io, "Minimal noise level=$(round(100 * mean(confidences .< 0.01), sigdigits=3))%. ",
            "Expected noise level=$(round(100 * mean(1 .- confidences), sigdigits=2))%.")
        println(io, "<br>")

        println(io, "<hr>\n<h1 id=\"Gene_Structure\">Gene structure</h1><br>")
        show(io, MIME("text/html"), n_tr_plot)
        println(io, "<br>")
        show(io, MIME("text/html"), gene_structure_plot)
        println(io, "</body>")
        println(io, "</html>")
    end

    @info "All done!"

    close(log_file)

    return 0
end

## All

run_cli(args::String) = run_cli(String.(Base.split(args)))

function run_cli(args::Vector{String}=ARGS)::Cint
    help_message = "Usage: baysor <command> [options]\n\nCommands:\n\trun\t\trun segmentation of the dataset\n\tpreview\t\tgenerate preview diagnostics of the dataset\n"

    debug = false
    if "--debug" in args
        args = args[args .!= "--debug"]
        debug = true
    end

    try
        if (length(args) == 0) || (length(args) == 1) && (args[1] == "-h" || args[1] == "--help")
            println(help_message)
            return 0
        end

        if args[1] == "run"
            return run_cli_main(args[2:end])
        end

        if args[1] == "preview"
            return run_cli_preview(args[2:end])
        end

        @error "Can't parse argument $(args[1])"
        println(help_message)
        return 1
    catch err
        if debug
            rethrow()
        else
            @error("$err\n\n" * join(["$s" for s in stacktrace(catch_backtrace())], "\n"))
        end
    end

    return 2
end

julia_main()::Cint = run_cli()