# Examples

Protocol-specific runnable examples live in the repository under
[examples/](https://github.com/kharchenkolab/Baysor/tree/cpp-0.9.0/examples):

- [Xenium pancreas](https://github.com/kharchenkolab/Baysor/tree/cpp-0.9.0/examples/Xenium_pancreas_membrane_377) —
  manifest input, molecule-label priors and Xenium Ranger import
- [ISS](https://github.com/kharchenkolab/Baysor/tree/cpp-0.9.0/examples/iss) —
  CSV molecule table with a prior workflow
- [osm-FISH](https://github.com/kharchenkolab/Baysor/tree/cpp-0.9.0/examples/osm-FISH) —
  image-mask prior workflow
- [STARmap](https://github.com/kharchenkolab/Baysor/tree/cpp-0.9.0/examples/STARmap) —
  3D dataset with per-layer polygon output

Each README covers data preparation and `baysor run` commands. The commands
use repository-relative config paths; run from the example directory in a
checkout, or download the matching [preset](configuration.md#protocol-presets)
and adjust the paths. For your own data, start with [Cell segmentation](run.md).
