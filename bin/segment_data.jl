using ArgParse
using Distributed
using Statistics

import CSV

import Baysor

function parse_commandline(args::Union{Nothing, Array}=nothing)
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
            default = "cell.assignment.csv"
        "--plot", "-p"
            help = "Save pdf with plot of the segmentation"
            action = :store_true

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
            help = "Standard deviation of center's prior distribution"
            arg_type = Float64
            default = 10.0
        "--shape-deg-freedom" # TODO: make it depend on mean number of molecules # TO JSON
            help = "Number of degrees of freedom for shape prior. Normally should be several times larger than expected number of molecules per cell."
            arg_type = Int
            default = 1000
        "--n-frames", "-n"
            help = "Number of frames, which is the same as number of processes. Algorithm data is splitted by frames to allow parallel run over frames."
            arg_type = Int
            default=1

        "--scale", "-s"
            help = "Scale parameter, which suggest approximate cell radius for the algorithm"
            arg_type = Float64
            required = true

        "coordinates"
            help = "CSV file with coordinates of transcripts and gene type"
            required = true
        "centers"
            help = "CSV file with coordinates of cell centers, extracted from DAPI staining"
            required = false
    end

    r = (args === nothing) ? parse_args(s) : parse_args(args, s)

    for k in ["gene", "x", "y"]
        r[k] = Symbol(r[k])
    end
    return r
end

load_df(args) = Baysor.load_df(args["coordinates"]; x_col=args["x"], y_col=args["y"], gene_col=args["gene"], min_molecules_per_gene=args["min-molecules-per-gene"])

function main()
    args = parse_commandline()

    @info "Run"
    @info "Load data..."
    df_spatial, gene_names = load_df(args);
    dfs_spatial = Baysor.split_spatial_data(df_spatial, args["n-frames"]);

    @info "Mean number of molecules per frame: $(median(size.(dfs_spatial, 1)))"

    @info "Done."

    size_prior = Baysor.ShapePrior(args["shape-deg-freedom"], [args["scale"], args["scale"]].^2);

    bm_data = nothing
    if args["centers"] !== nothing
        error("Not implemented yet")
        # df_centers = CSV.read(args["centers"]);
        # bm_data_prior = Baysor.initial_distributions(df_spatial, df_centers, args["center-std"]; size_prior=size_prior);
        # bm_data = Baysor.run_bmm(bm_data_prior, n_iters=args["iters"], new_component_weight=args["new-component-weight"],
        #                     prior_deg_freedom=args["shape-deg-freedom"]);
    else
        initial_params_per_frame = Baysor.cell_centers_with_clustering.(dfs_spatial, max(div(args["num-cells-init"], length(dfs_spatial)), 2); cov_mult=2);
        bm_data_arr = Baysor.initial_distributions.(dfs_spatial, initial_params_per_frame, size_prior=size_prior, new_component_weight=args["new-component-weight"]);

        addprocs(length(bm_data_arr) - 1)
        eval(:(@everywhere using Baysor))

        bm_data_res = Baysor.run_bmm_parallel(bm_data_arr, args["iters"], new_component_frac=args["new-component-fraction"],
                                              min_molecules_per_cell=args["min-molecules-per-cell"], n_refinement_iters=args["refinement-iters"]);

        bm_data = Baysor.merge_bm_data(bm_data_res);
    end

    if args["plot"]
        Baysor.plot_num_of_cells_per_iterarion(bm_data.tracer)

        # Plots.plot(
        #     Baysor.plot_num_of_cells_per_iterarion(bm_data.tracer),
        #     Baysor.plot_num_of_cells_per_iterarion(bm_data.tracer),
        #     layout=(2, 1)
        # )
        Plots.savefig("$(splitext(args["output"])[1]).pdf")
    end

    @info "Processing complete."

    @info "Save data to $(args["output"])"
    df_spatial[:assignment] = bm_data.assignment;
    df_spatial[:gene] = gene_names[df_spatial[:gene]]
    CSV.write(args["output"], df_spatial);
    @info "All done!"
end

main()