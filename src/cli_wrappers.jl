using ArgParse
using DataFrames
using Distributed
using ProgressMeter
using Statistics

import CSV

function parse_commandline(args::Union{Nothing, Array{String, 1}}=nothing)
    s = ArgParseSettings()
    @add_arg_table s begin
        "--x", "-x" # REPEAT IN JSON
            help = "name of x column"
            default = "x"
        "--y", "-y" # REPEAT IN JSON
            help = "name of gene column"
            default = "y"
        "--gene" # REPEAT IN JSON
            help = "name of gene column"
            default = "gene"

        "--iters", "-i"
            help = "Number of iterations"
            arg_type = Int
            default = 100
        "--refinement-iters" # TO JSON
            help = "Number of iterations for refinement of results"
            arg_type = Int
            default = 50
        "--min-molecules-per-gene" # TO JSON
            help = "Minimal number of molecules per gene"
            arg_type = Int
            default = 1
        "--min-molecules-per-cell" # TO JSON
            help = "Minimal number of molecules for a cell to be considered as real"
            arg_type = Int
            default = 3
        "--num-cells-init" # TO JSON
            help = "Initial number of cells. Ignored if CSV with centers is provided."
            arg_type = Int
            default = 100
        "--output", "-o"
            help = "Name of the output file"
            default = "segmentation.csv"
        "--plot", "-p"
            help = "Save pdf with plot of the segmentation"
            action = :store_true
        "--plot-frame-size"
            help = "Size of frame, which is used for result plotting. Ignored without '-v' option."
            arg_type = Float64
            default = 5000.0

        "--center-component-weight" # TO JSON
            help = "Prior weight of assignment a molecule to new component, created from DAPI centers. Paramter of Dirichlet process. Ignored if CSV with centers is not provided."
            arg_type = Float64
            default = 1.0
        "--new-component-weight" # TO JSON
            help = "Prior weight of assignment a molecule to new component. Paramter of Dirichlet process."
            arg_type = Float64
            default = 0.1
        "--new-component-fraction" # TO JSON
            help = "Fraction of distributions, sampled at each stage. Paramter of Dirichlet process."
            arg_type = Float64
            default = 0.2
        "--center-std" # TODO: remove it
            help = "Standard deviation of center's prior distribution. By default equal to scale."
            arg_type = Float64
        "--n-degrees-of-freedom-center" # TO JSON
            help = "Number of degrees of freedom for cell center distribution, used for posterior estimates of parameters. Ignored if centers are not provided. Default: equal to min-molecules-per-cell."
            arg_type = Int
        "--shape-deg-freedom" # TODO: make it depend on mean number of molecules # TO JSON
            help = "Number of degrees of freedom for shape prior. Normally should be several times larger than expected number of molecules per cell."
            arg_type = Int
            default = 1000
        "--n-frames", "-n"
            help = "Number of frames, which is the same as number of processes. Algorithm data is splitted by frames to allow parallel run over frames."
            arg_type = Int
            default=1
        "--gene-composition-neigborhood"
            help = "Number of neighbors (i.e. 'k'), which is used for gene composition visualization. Larger numbers leads to more global patterns."
            arg_type = Int
            default=20

        "--scale", "-s"
            help = "Scale parameter, which suggest approximate cell radius for the algorithm"
            arg_type = Float64

        "coordinates"
            help = "CSV file with coordinates of transcripts and gene type"
            required = true
        "centers"
            help = "CSV file with coordinates of cell centers, extracted from DAPI staining"
    end

    r = (args === nothing) ? parse_args(s) : parse_args(args, s)

    for k in ["gene", "x", "y"]
        r[k] = Symbol(r[k])
    end

    if r["center-std"] === nothing
        r["center-std"] = r["scale"]
    end

    if r["n-degrees-of-freedom-center"] === nothing
        r["n-degrees-of-freedom-center"] = r["min-molecules-per-cell"]
    end

    if r["centers"] === nothing && r["scale"] === nothing
        print("Either `centers` or `scale` must be provided.\n" * usage_string(s) * "\n")
        exit(1)
    end

    return r
end

load_df(args::Dict) = load_df(args["coordinates"]; x_col=args["x"], y_col=args["y"], gene_col=args["gene"], min_molecules_per_gene=args["min-molecules-per-gene"])

function plot_results(df_res::DataFrame, df_centers::Union{DataFrame, Nothing}, tracer::Dict, args::Dict{String,Any})
    ## Convergence
    p1 = plot_num_of_cells_per_iterarion(tracer);
    Plots.savefig("$(splitext(args["output"])[1])_convergence.pdf")

    ## Transcripts
    neighb_cm = neighborhood_count_matrix(df_res, args["gene-composition-neigborhood"]);
    color_transformation = gene_composition_transformation(neighb_cm)

    frame_size = args["plot-frame-size"]
    borders = [collect(range([floor(f(df_res[!, s]) / frame_size) * frame_size for f in [minimum, maximum]]..., step=frame_size)) for s in [:x, :y]];
    borders = collect(Iterators.product(borders...));

    plot_info = @showprogress "Extracting plot info..." pmap(borders) do b
        extract_plot_information(df_res, df_centers, b..., color_transformation=color_transformation, k=args["gene-composition-neigborhood"])
    end;
    plot_info = plot_info[length.(plot_info) .> 0];

    plot_width = 600
    p1 = Plots.plot([d[:plot] for d in plot_info]..., layout=(length(plot_info), 1), size=(plot_width, plot_width * length(plot_info)));
    Plots.savefig("$(splitext(args["output"])[1])_borders.pdf")
end

function run_cli(args::Union{Nothing, Array{String, 1}, String}=nothing)
    if args == "build"
        return 0
    end

    args = parse_commandline(args)

    @info "Run"
    @info "Load data..."
    df_spatial, gene_names = load_df(args)
    df_centers = args["centers"] === nothing ? nothing : read_spatial_df(args["centers"], x_col=args["x"], y_col=args["y"], gene_col=nothing)

    bm_data_arr = Baysor.initial_distribution_arr(df_spatial=df_spatial; n_frames=args["n-frames"],
        shape_deg_freedom=args["shape-deg-freedom"], scale=args["scale"], n_cells_init=args["num-cells-init"],
        new_component_weight=args["new-component-weight"], df_centers=df_centers, center_std=args["center-std"],
        center_component_weight=args["center-component-weight"], n_degrees_of_freedom_center=args["n-degrees-of-freedom-center"],
        min_molecules_per_cell=args["min-molecules-per-cell"]);

    addprocs(length(bm_data_arr) - nprocs())
    eval(:(@everywhere using Baysor))

    bm_data = run_bmm_parallel(bm_data_arr, args["iters"], new_component_frac=args["new-component-fraction"],
                               min_molecules_per_cell=args["min-molecules-per-cell"], n_refinement_iters=args["refinement-iters"]);

    @info "Processing complete."

    df_res = bm_data.x
    df_res[!, :assignment] = bm_data.assignment;

    if args["plot"]
        @info "Plot results"
        plot_results(df_res, df_centers, bm_data.tracer, args)
    end

    @info "Save data to $(args["output"])"
    df_res[!, :gene] = gene_names[df_res[!, :gene]]
    CSV.write(args["output"], df_res);

    @info "All done!"

    return 0
end