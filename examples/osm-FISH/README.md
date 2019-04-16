# Example run on osm-FISH data

Data is published in [S. Codeluppi, et al. - Spatial organization of the somatosensory cortex revealed by osmFISH](https://doi.org/10.1038/s41592-018-0175-z).

## Get the data

Download molecule coordinates in csv format (were extracted from the published [osmFISH raw coords](https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/osmFISH/data/mRNA_coords_raw_counting.hdf5) hdf5 file)

```bash
wget http://pklab.med.harvard.edu/viktor/spatial/baysor/mRNA_coords_raw_counting.csv
```

Download center coordinates in csv format (were extracted from the published [segmented regions](https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/osmFISH/data/polyT_seg.pkl) pkl file)

```bash
wget http://pklab.med.harvard.edu/viktor/spatial/baysor/centers_from_segmentation.csv
```

## CLI run

### Without DAPI

```bash
../../bin/segment_data.jl -i 500 -n 20 -p -s 50 ./mRNA_coords_raw_counting.csv
```

### With DAPI

```bash
../../bin/segment_data.jl -i 500 -n 20 -p -s 50 --center-component-weight 3.0 ./mRNA_coords_raw_counting.csv ./centers_from_segmentation.csv
```