# Segmentation

## Normal run

To run the algorithm on your data, use the following command:

```bash
baysor run <args> [options] [flags]
```

CLI parameters:

```@docs
Baysor.CLI.run
```

For the description of all config parameters, see [example_config.toml](https://github.com/kharchenkolab/Baysor/blob/master/configs/example_config.toml).

## Using a prior segmentation

In some cases, you may want to use another segmentation as a prior for Baysor. The most popular case is having a segmentation based on DAPI/poly-A stainings: such information helps to understand where nuclei are positioned, but it's often quite imprecise. To take this segmentation into account you can pass it as the second positional argument to Baysor:

```bash
baysor run [ARGS] MOLECULES_FILE [PRIOR_SEGMENTATION]
```

Here, `PRIOR_SEGMENTATION` can be a path to a binary image with a segmentation mask, an image with integer cell segmentation labels or a column name in the `MOLECULES_FILE` with integer cell assignment per molecule (`0` value means no assignment). In the latter case, the column name must have `:` prefix, e.g. for column `cell` you should use `baysor run [ARGS] molecules.csv :cell`. In case the image is too big to be stored in the tiff format, Baysor supports MATLAB '.mat' format: it should contain a single field with an integer matrix for either a binary mask or segmentation labels. When loading the segmentation, Baysor filters segments that have less than `min-molecules-per-segment` molecules. It can be set in the toml config, and the default value is `min-molecules-per-segment = min-molecules-per-cell / 4`. **Note:** only CSV column prior is currently supported for 3D segmentation.

To specify the expected quality of the prior segmentation you may use `prior-segmentation-confidence` parameter. The value `0.0` makes the algorithm ignore the prior, while the value `1.0` restricts the algorithm from contradicting the prior. Prior segmentation is mainly needed for the cases where gene expression signal is not enough, e.g. with very sparse protocols (such as ISS or DARTFISH). Another potential use case is high-quality data with a visible sub-cellular structure. In these situations, setting `prior-segmentation-confidence > 0.7` is recommended. Otherwise, the default value `0.2` should work well.

### Segmenting stains

If you have a non-segmented DAPI image, the simplest way to segment it would go through the following steps [ImageJ](https://fiji.sc/):

1. Open the image (File -> Open)
2. Go to Image -> Type and pick "8-bit"
3. Run Process -> Filters -> Gaussian Blur, using Sigma = 1.0. *The value can vary, depending on your DAPI, but 1.0 generally works fine.*
4. Run Image -> Adjust -> Auto Threshold, using Method = Default. *Different methods can give the best results for different cases. Often "Mean" also works well.*
5. Run Process -> Binary -> Watershed
6. Save the resulting image in the .tif

Another promising tool is [CellPose](https://github.com/mouseland/cellpose), however, it may require some manual labeling to fine-tune the network.

## Segmenting cells with pronounced intracellular structure

High-resolution protocols, such as MERFISH or seq-FISH, can capture the intracellular structure. Most often, it would mean a pronounced difference between nuclear and cytoplasmic gene composition. By default, such differences would push Baysor to recognize compartments as different cells. However, if some compartment-specific genes are known, they may be used to mitigate the situation. These genes can be specified through `--config.segmentation.nuclei-genes` and `--config.segmentation.cyto-genes` options, *e.g.*:

```julia
baysor run -m 30 --n-clusters=1 -s 30 --scale-std=50% --config.segmentation.nuclei-genes=Neat1 --config.segmentation.cyto-genes=Apob,Net1,Slc5a1,Mptx2 --config.data.exclude-genes='Blank*' ./molecules.csv
```

Please, notice that it's highly recommended to set `--n-clusters=1`, so molecule clustering would not be affected by compartment differences.

**Note.** Currently, there is no automated way to determine such compartment-specific genes. So, the only way we can suggest is interactive explaration of data. In theory, it should be straightforward to infer such information from DAPI and poly-A stains, however, it is not implemented yet. If you have a particular need for such functionality, please submit an issue with the description of your experimental setup.

## Outputs

### Segmentation results

- **segmentation\_counts.loom** or **segmentation\_counts.tsv** (depends on `--count-matrix-format`): count matrix with segmented stats. In the case of loom format, column attributes also contain the same info as **segmentation\_cell\_stats.csv**.
- **segmentation.csv**: segmentation info per molecule:
  - `confidence`: probability of a molecule to be real (i.e. not noise)
  - `cell`: id of the assigned cell. Value "" corresponds to noise.
  - `cluster`: id of the molecule cluster
  - `assignment_confidence`: confidence that the molecule is assigned to a correct cell
  - `is_noise`: shows whether molecule was assigned to noise *(it equals `true` if and only if `cell` == "")*
  - `ncv_color`: RGB code of the neighborhood composition coloring
- **segmentation\_cell\_stats.csv**: diagnostic info about cells. The following parameters can be used to filter low-quality cells:
  - `area`: area of the convex hull around the cell molecules
  - `avg_confidence`: average confidence of the cell molecules
  - `density`: the number of molecules in a cell divided by the cell area
  - `elongation`: ratio of the two eigenvalues of the cell covariance matrix
  - `n_transcripts`: number of molecules per cell
  - `avg_assignment_confidence`: average assignment confidence per cell. Cells with low `avg_assignment_confidence` have a much higher chance of being an artifact.
  - `max_cluster_frac` *(only if `n-clusters > 1`)*: fraction of the molecules coming from the most popular cluster. Cells with low `max_cluster_frac` are often doublets.
  - `lifespan`: number of iterations the given component exists. The maximal `lifespan` is clipped proportionally to the total number of iterations. Components with a short lifespan likely correspond to noise.

### Visualization

- **segmentation\_polygons\_2d/3d.json**: polygons used for visualization in GeoJSON format. In the case of 3D segmentation, `2d.json` file contains polygons for all molecules pulled across the z-stack. And 3D shows polygons per z-slice. In case of continuous z, it's binned into 20 uniform bins. Depending on the format, it can contain `GeometryCollection` or `FeatureCollection`:
    - For 3D, the file contains an array of dictionaries (one per z-slice), each of which representing a `Collection`. For 2D data it's just a single dictionary with a `Collection`.
    - Each `GeometryCollection` has a field `geometries`, which is an array of polygons with `cell` field set to cell ids and `coordinates` set to its coordinates.
    - `FeatureCollection` is the format, compatible with [10x SpaceRanger](https://www.10xgenomics.com/support/software/xenium-ranger/1.7/analysis/inputs/XR-input-overview#compat-files). It contains a list of `Feature`s with cell ids saved in the `id` field and coordinates in `geometry/coordinates`.
- **segmentation\_diagnostics.html**: visualization of the algorithm QC. *Shown only when `-p` is set.*
- **segmentation\_borders.html**: visualization of cell borders for the dataset colored by local gene expression composition (first part) and molecule clusters (second part). *Shown only when `-p` is set.*

### Other

- **segmentation\_params.dump.toml**: aggregated parameters from the config and CLI

## Choice of parameters

Most important parameters:

- `scale` is the most sensitive parameter, which specifies the expected radius of a cell. It doesn't have to be precise, but the wrong setup can lead to over- or under-segmentation. This parameter is inferred automatically if cell centers are provided.
- `min-molecules-per-cell` is the number of molecules, required for a cell to be considered as real. It really depends on the protocol. For instance, for ISS it's fine to set it to 3, while for MERFISH it can require hundreds of molecules.

Run parameters:

- `--config.segmentation.n-cells-init` expected number of cells in data. This parameter influences the convergence speed of the algorithm, as well as peak memory usage. Setting this value too small would lead to under-segmentation.
- `--config.segmentation.iters` number of iterations for the algorithm. **At the moment, no convergence criteria are implemented, so it will work exactly `iters` iterations**. Thus, too small values would lead to non-convergence of the algorithm, while larger ones would just increase working time. Optimal values can be estimated by the convergence plots in the diagnostics report.
