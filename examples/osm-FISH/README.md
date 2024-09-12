# Example run on osm-FISH data

Data is published in [S. Codeluppi, et al. - Spatial organization of the somatosensory cortex revealed by osmFISH](https://doi.org/10.1038/s41592-018-0175-z).

## Get the data

Create directories:

```bash
mkdir -p data output_dapi output_no_dapi
```

Download molecule coordinates in the csv format (were extracted from the published [osmFISH raw coords](https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/osmFISH/data/mRNA_coords_raw_counting.hdf5) hdf5 file)

```bash
wget -P data http://pklab.med.harvard.edu/viktor/baysor/osm_fish/mRNA_coords_raw_counting.csv
```

Download center coordinates in csv format (were extracted from the published [segmented regions](https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/osmFISH/data/polyT_seg.pkl) pkl file)

```bash
wget -P data http://pklab.med.harvard.edu/viktor/baysor/osm_fish/centers_from_segmentation.csv
```

Or convert it to an image using the following python code:

```python
from skimage import io
import numpy as np
import pandas as pd

segmentation_labels = pd.read_pickle("./data/polyT_seg.pkl")
df = pd.read_csv("./data/mRNA_coords_raw_counting.csv")

dapi_shape = np.concatenate([[x.max(axis=0)] for x in segmentation_labels.values()]).max(axis=0)
dapi_shape = np.maximum(dapi_shape, df[["x", "y"]].max().astype(int).values)
segmentation_mask = np.zeros(dapi_shape + 1, dtype=int)

for label,coords in segmentation_labels.items():
    segmentation_mask[coords[:,0], coords[:,1]] = int(label)

io.imsave("./data/segmentation.tiff", segmentation_mask.T.astype(np.uint16))
```

## CLI run

WARNING: these command run the segmentation using 20 threads. To reduce this number change `JULIA_NUM_THREADS` option.

### Without DAPI

```bash
JULIA_NUM_THREADS=20 baysor -i 500 -c ../../configs/osm_fish.toml -p -o ./output_no_dapi ./data/mRNA_coords_raw_counting.csv
```

### With DAPI

```bash
JULIA_NUM_THREADS=20 baysor  -i 500 -p -c ../../configs/osm_fish.toml -o ./output_dapi -p ./data/mRNA_coords_raw_counting.csv ./data/centers_from_segmentation.csv
```