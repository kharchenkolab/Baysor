using Configurations

@option mutable struct DataOptions # --config.data.*
    x::String = "x" # Name of the x column in the input data. Default: "x"
    y::String = "y" # Name of the y column in the input data. Default: "y"
    z::String = "z" # Name of the y column in the input data. Default: "z"
    gene::String = "gene" # Name of gene column in the input data. Default: "gene"
    force_2d::Bool = false # Ignores z-column in the data if it is provided
    min_molecules_per_gene::Int = 1 # Minimal number of molecules per gene. Default: 1
    exclude_genes::String = "" # Comma-separated list of genes or regular expressions to ignore during segmentation. Example: 'Blank*,MALAT1'
    min_molecules_per_cell::Int = 0 # Minimal number of molecules for a cell to be considered as real. It's an important parameter, as it's used to infer several other parameters. Default: 3
    min_molecules_per_segment::Int = 0 # Minimal number of molecules in a segmented region, required for this region to be considered as a possible cell. Default: min-molecules-per-cell / 4
    confidence_nn_id = 0 # Number of nearest neighbors to use for confidence estimation. Default: min-molecules-per-cell / 2 + 1
end

@option mutable struct SegmentationOptions # --config.segmentation.*
    scale::Float64 = -1.0 # Negative values mean it must be estimated from `min_molecules_per_cell`
    scale_std::String = "25%" # Standard deviation of scale across cells. Can be either number, which means absolute value of the std, or string ended with "%" to set it relative to scale. Default: "25%"
    estimate_scale_from_centers::Bool = true # Use scale estimate from DAPI if provided. Default: true

    n_clusters::Int = 4 # Number of clusters to use for cell type segmentation. Default: 4
    prior_segmentation_confidence::Float64 = 0.2 # Confidence of the prior segmentation. Default: 0.2
    iters::Int = 500 # Number of iterations for the cell segmentation algorithm. Default: 500
    n_cells_init::Int = 0 # Initial number of cells

    nuclei_genes::String = "" # Comma-separated list of nuclei-specific genes. If provided, `cyto-genes` has to be set, as well.
    cyto_genes::String = "" # Comma-separated list of cytoplasm-specific genes. If provided, `nuclei-genes` has to be set, as well.

    # The parameters below are not supposed to be changed normally
    new_component_weight::Float64 = 0.2 # Prior weight of assignment a molecule to new component. Default: 0.2
    new_component_fraction::Float64 = 0.3 # Fraction of distributions, sampled at each stage. Default: 0.3
end

@option mutable struct PlottingOptions # --config.plotting.*
    gene_composition_neigborhood::Int = 0 # Number of neighbors (i.e. 'k' in k-NN), which is used for gene composition visualization. Larger numbers leads to more global patterns. Default: estimate from min-molecules-per-cell
    min_pixels_per_cell::Int = 15 # Number of pixels per cell of minimal size, used to estimate size of the final plot. For most protocols values around 7-30 give enough visualization quality. Default: 15
end

@option mutable struct RunOptions # --config.*
    data::DataOptions = DataOptions()
    segmentation::SegmentationOptions = SegmentationOptions()
    plotting::PlottingOptions = PlottingOptions()
end