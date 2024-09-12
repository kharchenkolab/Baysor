# Example run on pciSeq in-situ sequencing (ISS) data

Data is published in [Qian X, et al. - A spatial atlas of inhibitory cell types in mouse hippocampus](https://doi.org/10.1101/431957) and available at [figshare.com](https://figshare.com/s/88a0fc8157aca0c6f0e8).

## Get the data

Create directories:

```bash
mkdir -p data output_dapi output_no_dapi
```

And either download raw data and pre-process it, or download already pre-processed data in the compatible format.

### Preprocessing raw data

In this example we used "3-3" section. So, you need to download "pciSeq_3-3_right.mat" file ("pciSeq_3-3_left.mat" works as well, because they have the same set of molecules) and the corresponding DAPI ("DAPI_3-3.jpg"). To convert .mat file to the compatible csv format, you need to define ISS class (download [the file](https://github.com/kdharris101/iss/blob/a2fab5ec452d447f5b8adfde5a3d85be49aa25be/%40iss/iss.m) and run `run("iss.m")` inside your MATLAB session. After that, run the following code:

```matlab
load("pciSeq_3-3_right.mat")

y = o.SpotGlobalYX(:, 1);
x = o.SpotGlobalYX(:, 2);
gene = o.GeneNames(o.SpotCodeNo);
is_combinatorial = o.SpotCombi;
intensity = o.SpotIntensity;
score = o.SpotScore;

T = table(x, y, gene, is_combinatorial, intensity, score);
writetable(T, "pciSeq_3-3.csv");
```

<!--
[P,I] = max(o.pSpotCell, [], 2) - put cell assignment to I
o.CellYX - cell centers
But these are subset of CA1 section (right and left correspondingly).
Spots, stored in left and right files are identical and global for all section.
 -->

Afterwards, to obtain segmentation from the DAPI, we ran Watershed segmentation using ImageJ:

1. Open DAPI_3-3.jpg image
2. Run Process -> Filters -> Gaussian Blur, using Sigma = 1.0
3. Run Image -> Adjust -> Auto Threshold, using Method = Mean
4. Run Process -> Binary -> Watershed
5. Save the image in the .tif format using name "DAPI_3-3_mask.tif"

### Download pre-processed data

```bash
wget -P data http://pklab.med.harvard.edu/viktor/baysor/iss/DAPI_3-3_mask.tif
wget -P data http://pklab.med.harvard.edu/viktor/baysor/iss/pciSeq_3-3.csv
```

## CLI run

WARNING: these command run the segmentation using 20 threads. To reduce this number change `JULIA_NUM_THREADS` option.

### Without DAPI

```bash
JULIA_NUM_THREADS=20 baysor -i 500 -p -c ../../configs/iss.toml -o ./output_no_dapi -p ./data/pciSeq_3-3.csv
```

### With DAPI

```bash
JULIA_NUM_THREADS=20 baysor -i 500 -p -c ../../configs/iss.toml -o ./output_dapi -p ./data/pciSeq_3-3.csv ./data/DAPI_3-3_mask.tif
```