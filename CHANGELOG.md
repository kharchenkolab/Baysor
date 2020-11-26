## Upcoming

### Fixed

- Small bug fixes
- `min-molecules-per-segment` parameter is working now
- Fixed visualization of prior segmentation

## [0.4.2] — 2020-11-26

### Changed

- Fixed Makefile julia version
- Improved polygon visualization
- Regressed to Plots 1.6.0 because of the performance issues
- Fixed docker build

## [0.4.1] — 2020-10-30

### Added

- Added travis config
- Saving NCV colors to the `ncv_color` field of *segmentation.csv*

### Removed

- Dropped support for julia < 1.5

### Changed

- Updated dependencies
- `find_grid_point_labels_kde` now preserves label ids
- Fixed docker build
- Added ImageMagick dependency to fix problems with DAPI prior
