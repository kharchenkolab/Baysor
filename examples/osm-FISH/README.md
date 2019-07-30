# Example run on osm-FISH data

Data is published in [S. Codeluppi, et al. - Spatial organization of the somatosensory cortex revealed by osmFISH](https://doi.org/10.1038/s41592-018-0175-z).

## Get the data

Create directories:

```bash
mkdir -p data output_dapi output_no_dapi
```

Download molecule coordinates in the csv format (were extracted from the published [osmFISH raw coords](https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/osmFISH/data/mRNA_coords_raw_counting.hdf5) hdf5 file)

```bash
wget -P data http://pklab.med.harvard.edu/viktor/spatial/baysor/mRNA_coords_raw_counting.csv
```

Download center coordinates in csv format (were extracted from the published [segmented regions](https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/osmFISH/data/polyT_seg.pkl) pkl file)

```bash
wget -P data http://pklab.med.harvard.edu/viktor/spatial/baysor/centers_from_segmentation.csv
```

## CLI run

WARNING: these command run the segmentation using 20 proccesses. To reduce this number change `-n 20` option.

### Without DAPI

```bash
../../bin/segment_data.jl -i 500 --num-cells-init=30000 -n 20 -c ../../configs/osm_fish.toml -p -o ./output_no_dapi ./data/mRNA_coords_raw_counting.csv
```

### With DAPI

```bash
../../bin/segment_data.jl  -i 500 --num-cells-init=30000 -n 20 -p -c ../../configs/osm_fish.toml -o ./output_dapi -p ./data/mRNA_coords_raw_counting.csv ./data/centers_from_segmentation.csv
```