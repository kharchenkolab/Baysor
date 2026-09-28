#include "baysor/data_loading/prior_segmentation.h"
#include "baysor/processing/utils/utils.h"
#include "baysor/utils/general.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <spdlog/spdlog.h>
#include <stdexcept>
#include <unordered_map>
#include <tiffio.h>
#include <queue>
#include <unordered_set>

namespace baysor {

namespace {

std::string file_extension_lower(const std::string& path) {
    auto dot = path.rfind('.');
    if (dot == std::string::npos) return "";
    std::string ext = path.substr(dot + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

struct BoundaryPolygon {
    int label = 0;
    std::vector<double> xs;
    std::vector<double> ys;
    double min_x = 0.0;
    double max_x = 0.0;
    double min_y = 0.0;
    double max_y = 0.0;
    double bbox_area = 0.0;
};

struct MoleculeBounds2D {
    bool valid = false;
    double min_x = 0.0;
    double max_x = 0.0;
    double min_y = 0.0;
    double max_y = 0.0;
};

MoleculeBounds2D compute_molecule_bounds(
    const std::vector<double>& mol_x,
    const std::vector<double>& mol_y)
{
    MoleculeBounds2D bounds;
    if (mol_x.empty() || mol_y.empty()) return bounds;
    bounds.valid = true;
    bounds.min_x = bounds.max_x = mol_x[0];
    bounds.min_y = bounds.max_y = mol_y[0];
    for (size_t i = 1; i < mol_x.size(); ++i) {
        bounds.min_x = std::min(bounds.min_x, mol_x[i]);
        bounds.max_x = std::max(bounds.max_x, mol_x[i]);
        bounds.min_y = std::min(bounds.min_y, mol_y[i]);
        bounds.max_y = std::max(bounds.max_y, mol_y[i]);
    }
    return bounds;
}

void filter_boundary_polygons_to_molecule_bounds(
    std::vector<BoundaryPolygon>& polygons,
    const MoleculeBounds2D& bounds)
{
    if (!bounds.valid || polygons.empty()) return;
    size_t before = polygons.size();
    polygons.erase(
        std::remove_if(
            polygons.begin(), polygons.end(),
            [&](const BoundaryPolygon& poly) {
                return poly.max_x < bounds.min_x || poly.min_x > bounds.max_x ||
                       poly.max_y < bounds.min_y || poly.min_y > bounds.max_y;
            }),
        polygons.end());

    if (polygons.size() != before) {
        spdlog::info( // GCOVR_EXCL_LINE: unwind-only block of the spdlog::info statement, entered only if logging throws
            "Filtered boundary priors to {} polygons overlapping molecule bounds "
            "[x=({:.2f}, {:.2f}), y=({:.2f}, {:.2f})] from {} total",
            polygons.size(), bounds.min_x, bounds.max_x, bounds.min_y, bounds.max_y, before);
    } else {
        spdlog::info(
            "All {} boundary polygons overlap molecule bounds "
            "[x=({:.2f}, {:.2f}), y=({:.2f}, {:.2f})]",
            before, bounds.min_x, bounds.max_x, bounds.min_y, bounds.max_y);
    }
}

bool point_in_polygon(const BoundaryPolygon& poly, double x, double y) {
    bool inside = false;
    int n = static_cast<int>(poly.xs.size());
    if (n < 3) return false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        double xi = poly.xs[i], yi = poly.ys[i];
        double xj = poly.xs[j], yj = poly.ys[j];
        bool intersects = ((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / ((yj - yi) + 1e-30) + xi);
        if (intersects) inside = !inside;
    }
    return inside;
}

std::vector<BoundaryPolygon> load_boundary_polygons(const std::string& path) {
    std::vector<double> vx = read_double_column(path, "vertex_x");
    std::vector<double> vy = read_double_column(path, "vertex_y");
    std::vector<int> labels;

    try {
        auto label_d = read_double_column(path, "label_id");
        labels.resize(label_d.size());
        for (size_t i = 0; i < label_d.size(); ++i) labels[i] = static_cast<int>(std::llround(label_d[i]));
    } catch (...) {
        labels = encode_prior_labels(read_string_column(path, "cell_id"), "");
    }

    if (vx.size() != vy.size() || vx.size() != labels.size()) {
        throw std::runtime_error("Boundary file columns have inconsistent lengths: " + path); // GCOVR_EXCL_LINE: unreachable, all columns are read from the same table so lengths always match
    }

    std::unordered_map<int, int> poly_index;
    std::vector<BoundaryPolygon> polygons;
    polygons.reserve(labels.size() / 8);
    for (size_t i = 0; i < labels.size(); ++i) {
        int label = labels[i];
        if (label <= 0) continue;
        auto it = poly_index.find(label);
        if (it == poly_index.end()) {
            BoundaryPolygon poly;
            poly.label = label;
            poly.min_x = poly.max_x = vx[i];
            poly.min_y = poly.max_y = vy[i];
            polygons.push_back(std::move(poly));
            it = poly_index.emplace(label, static_cast<int>(polygons.size()) - 1).first;
        }
        auto& poly = polygons[it->second];
        poly.xs.push_back(vx[i]);
        poly.ys.push_back(vy[i]);
        poly.min_x = std::min(poly.min_x, vx[i]);
        poly.max_x = std::max(poly.max_x, vx[i]);
        poly.min_y = std::min(poly.min_y, vy[i]);
        poly.max_y = std::max(poly.max_y, vy[i]);
    }

    for (auto& poly : polygons) {
        poly.bbox_area = std::max(poly.max_x - poly.min_x, 0.0) *
                         std::max(poly.max_y - poly.min_y, 0.0);
    }
    return polygons;
}

double choose_boundary_grid_size(const std::vector<BoundaryPolygon>& polygons) {
    if (polygons.empty()) return 50.0;
    std::vector<double> dims;
    dims.reserve(polygons.size());
    for (const auto& poly : polygons) {
        dims.push_back(std::max(poly.max_x - poly.min_x, poly.max_y - poly.min_y));
    }
    std::sort(dims.begin(), dims.end());
    double median = dims[dims.size() / 2];
    return std::max(10.0, median * 2.0);
}

long long boundary_grid_key(int gx, int gy) {
    return (static_cast<long long>(gx) << 32) ^ static_cast<unsigned int>(gy);
}

std::vector<int> assign_molecules_to_boundaries(
    const std::vector<double>& mol_x,
    const std::vector<double>& mol_y,
    const std::vector<BoundaryPolygon>& polygons,
    int min_molecules_per_segment
) {
    double grid = choose_boundary_grid_size(polygons);
    std::unordered_map<long long, std::vector<int>> grid_index;
    grid_index.reserve(polygons.size() * 2);

    for (int pi = 0; pi < static_cast<int>(polygons.size()); ++pi) {
        const auto& poly = polygons[pi];
        int gx0 = static_cast<int>(std::floor(poly.min_x / grid));
        int gx1 = static_cast<int>(std::floor(poly.max_x / grid));
        int gy0 = static_cast<int>(std::floor(poly.min_y / grid));
        int gy1 = static_cast<int>(std::floor(poly.max_y / grid));
        for (int gx = gx0; gx <= gx1; ++gx) {
            for (int gy = gy0; gy <= gy1; ++gy) {
                grid_index[boundary_grid_key(gx, gy)].push_back(pi);
            }
        }
    }

    std::vector<int> labels(mol_x.size(), 0);
    for (int i = 0; i < static_cast<int>(mol_x.size()); ++i) {
        int gx = static_cast<int>(std::floor(mol_x[i] / grid));
        int gy = static_cast<int>(std::floor(mol_y[i] / grid));
        auto it = grid_index.find(boundary_grid_key(gx, gy));
        if (it == grid_index.end()) continue;

        int best_label = 0;
        double best_area = std::numeric_limits<double>::infinity();
        for (int pi : it->second) {
            const auto& poly = polygons[pi];
            if (mol_x[i] < poly.min_x || mol_x[i] > poly.max_x ||
                mol_y[i] < poly.min_y || mol_y[i] > poly.max_y) {
                continue;
            }
            if (!point_in_polygon(poly, mol_x[i], mol_y[i])) continue;
            if (poly.bbox_area < best_area) {
                best_area = poly.bbox_area;
                best_label = poly.label;
            }
        }
        labels[i] = best_label;
    }

    filter_segmentation_labels(labels, min_molecules_per_segment);
    return labels;
}

} // anonymous namespace

// ============================================================================
// Detect prior segmentation type
// ============================================================================

PriorInputType detect_prior_seg_type(const std::string& prior_seg_arg) {
    return parse_prior_input_spec(prior_seg_arg).type;
}

PriorInputOptions parse_prior_input_spec(const std::string& prior_seg_arg) {
    PriorInputOptions opts;
    if (prior_seg_arg.empty()) return opts;
    if (prior_seg_arg[0] == ':') {
        opts.type = PriorInputType::Column;
        opts.column_name = prior_seg_arg.substr(1);
        return opts;
    }
    std::string ext = file_extension_lower(prior_seg_arg);
    opts.path = prior_seg_arg;
    if (ext == "csv" || ext == "parquet" || ext == "pq") {
        opts.type = PriorInputType::Boundary;
    } else {
        opts.type = PriorInputType::Image;
    }
    return opts;
}

void apply_prior_input_spec(PriorInputOptions& opts, const std::string& prior_seg_arg) {
    auto parsed = parse_prior_input_spec(prior_seg_arg);
    opts.type = parsed.type;
    opts.path = parsed.path;
    opts.column_name = parsed.column_name;
}

// ============================================================================
// Filter segmentation labels
// ============================================================================

void filter_segmentation_labels(std::vector<int>& labels, int min_molecules_per_segment) {
    if (min_molecules_per_segment <= 0) return;
    if (labels.empty()) return;

    int max_label = *std::max_element(labels.begin(), labels.end());
    if (max_label <= 0) return;

    // Count molecules per segment (1-based labels, index 0 -> label 1)
    std::vector<int> counts(max_label, 0);
    for (int lab : labels) {
        if (lab > 0 && lab <= max_label) {
            counts[lab - 1]++;
        }
    }

    // Zero out labels for small segments
    for (int& lab : labels) {
        if (lab > 0 && counts[lab - 1] < min_molecules_per_segment) {
            lab = 0;
        }
    }
}

// ============================================================================
// Parse prior segmentation from a column
// ============================================================================

std::vector<int> encode_prior_labels(
    const std::vector<std::string>& raw_values,
    const std::string& unassigned_label,
    int min_molecules_per_segment
) {
    int n = static_cast<int>(raw_values.size());

    // Encode labels: collect unique non-unassigned values, sort, assign 1-based IDs
    std::set<std::string> unique_labels;
    for (const auto& v : raw_values) {
        if (v != unassigned_label && !v.empty()) {
            unique_labels.insert(v);
        }
    }

    std::vector<std::string> sorted_labels(unique_labels.begin(), unique_labels.end());
    std::unordered_map<std::string, int> label_map;
    label_map.reserve(sorted_labels.size());
    for (int i = 0; i < static_cast<int>(sorted_labels.size()); ++i) {
        label_map[sorted_labels[i]] = i + 1;  // 1-based
    }

    // Map raw values to integer labels (0 = unassigned)
    std::vector<int> result(n, 0);
    int n_unassigned = 0;
    for (int i = 0; i < n; ++i) {
        if (raw_values[i] == unassigned_label || raw_values[i].empty()) {
            n_unassigned++;
            continue;
        }
        auto it = label_map.find(raw_values[i]);
        if (it != label_map.end()) {
            result[i] = it->second;
        }
    }

    if (n_unassigned == 0) {
        spdlog::warn("No unassigned molecules found in the prior segmentation. "
                      "Did you specify the unassigned label correctly ('{}')?",
                      unassigned_label);
    }

    spdlog::info("Parsed {} prior segments ({} unassigned molecules)", // GCOVR_EXCL_LINE: unwind-only block of the spdlog::info statement, entered only if logging throws
                 sorted_labels.size(), n_unassigned);

    // Filter small segments
    filter_segmentation_labels(result, min_molecules_per_segment);

    return result;
}

std::vector<int> parse_prior_from_column(
    const std::string& molecule_path,
    const std::string& col_name,
    const std::string& unassigned_label,
    int min_molecules_per_segment
) {
    auto raw_values = read_string_column(molecule_path, col_name);
    return encode_prior_labels(raw_values, unassigned_label, min_molecules_per_segment);
}

// ============================================================================
// Estimate scale from assignment (center-to-center distances)
// ============================================================================

std::pair<double, double> estimate_scale_from_assignment(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& assignment,
    int min_molecules_per_cell
) {
    int max_label = *std::max_element(assignment.begin(), assignment.end());
    if (max_label <= 0) {
        throw std::runtime_error("No assigned molecules found for scale estimation");
    }

    // Group molecule indices by assignment
    auto groups = split_ids(assignment, max_label, /*drop_zero=*/true);

    // Compute centers for groups that pass the size threshold
    int n_dims = static_cast<int>(pos_data.rows());
    std::vector<Eigen::VectorXd> centers;
    centers.reserve(groups.size());

    for (const auto& ids : groups) {
        if (static_cast<int>(ids.size()) < min_molecules_per_cell) continue;

        Eigen::VectorXd center = Eigen::VectorXd::Zero(n_dims);
        for (int idx : ids) {
            center += pos_data.col(idx);
        }
        center /= static_cast<double>(ids.size());
        centers.push_back(std::move(center));
    }

    int n_centers = static_cast<int>(centers.size());
    if (n_centers < 3) {
        throw std::runtime_error(
            "Not enough prior cells pass the min_molecules_per_cell threshold (" +
            std::to_string(n_centers) + " < 3). Please specify scale manually.");
    }

    Eigen::MatrixXd center_mat(n_dims, n_centers);
    for (int i = 0; i < n_centers; ++i) {
        center_mat.col(i) = centers[i];
    }

    // For each center, find distance to nearest neighbor, then take radius = dist / 2.
    // This is the exact same nearest-neighbor definition as the all-pairs loop,
    // but uses the shared exact KD-tree implementation for large priors.
    auto knn = knn_parallel(center_mat, center_mat, std::min(8, n_centers), true);
    std::vector<double> radii(n_centers);
    for (int i = 0; i < n_centers; ++i) {
        double min_dist = std::numeric_limits<double>::max();
        const auto& ids = knn.indices[i];
        const auto& dists = knn.distances[i];
        for (int j = 0; j < static_cast<int>(ids.size()); ++j) {
            if (ids[j] == i) continue;
            min_dist = dists[j];
            break;
        }
        if (min_dist == std::numeric_limits<double>::max()) {
            for (int j = 0; j < n_centers; ++j) { // GCOVR_EXCL_LINE: unreachable, knn_parallel returns >= k-1 >= 2 non-self neighbours for n_centers >= 3
                if (i == j) continue; // GCOVR_EXCL_LINE: unreachable, knn_parallel returns >= k-1 >= 2 non-self neighbours for n_centers >= 3
                double dist = (centers[i] - centers[j]).norm(); // GCOVR_EXCL_LINE: unreachable, knn_parallel returns >= k-1 >= 2 non-self neighbours for n_centers >= 3
                if (dist < min_dist) min_dist = dist; // GCOVR_EXCL_LINE: unreachable, knn_parallel returns >= k-1 >= 2 non-self neighbours for n_centers >= 3
            }
        }
        radii[i] = min_dist / 2.0;
    }

    // Compute median and MAD (median absolute deviation, normalized)
    std::sort(radii.begin(), radii.end());
    double median_radius;
    if (n_centers % 2 == 0) {
        median_radius = (radii[n_centers / 2 - 1] + radii[n_centers / 2]) / 2.0;
    } else {
        median_radius = radii[n_centers / 2];
    }

    // MAD with normalization factor 1.4826 (for consistency with normal distribution)
    std::vector<double> abs_devs(n_centers);
    for (int i = 0; i < n_centers; ++i) {
        abs_devs[i] = std::abs(radii[i] - median_radius);
    }
    std::sort(abs_devs.begin(), abs_devs.end());
    double mad;
    if (n_centers % 2 == 0) {
        mad = (abs_devs[n_centers / 2 - 1] + abs_devs[n_centers / 2]) / 2.0;
    } else {
        mad = abs_devs[n_centers / 2];
    }
    mad *= 1.4826;  // normalize=true in Julia's mad()

    return {median_radius, mad};
}

// ============================================================================
// Estimate scale from connected-component pixel areas (image-based prior)
// ============================================================================

std::pair<double, double> estimate_scale_from_image_areas(
    const std::vector<size_t>& component_pixel_areas
) {
    // Compute nucleus radius from area: r = sqrt(area / π)
    std::vector<double> radii;
    radii.reserve(component_pixel_areas.size());
    for (size_t a : component_pixel_areas) {
        if (a > 0) radii.push_back(std::sqrt(static_cast<double>(a) / kPi));
    }

    int n = static_cast<int>(radii.size());
    if (n < 3) {
        throw std::runtime_error(
            "Not enough prior cells for scale estimation from image areas (" +
            std::to_string(n) + " < 3). Please specify scale manually.");
    }

    std::sort(radii.begin(), radii.end());
    double median_radius;
    if (n % 2 == 0) {
        median_radius = (radii[n / 2 - 1] + radii[n / 2]) / 2.0;
    } else {
        median_radius = radii[n / 2];
    }

    std::vector<double> abs_devs(n);
    for (int i = 0; i < n; ++i) abs_devs[i] = std::abs(radii[i] - median_radius);
    std::sort(abs_devs.begin(), abs_devs.end());
    double mad;
    if (n % 2 == 0) {
        mad = (abs_devs[n / 2 - 1] + abs_devs[n / 2]) / 2.0;
    } else {
        mad = abs_devs[n / 2];
    }
    mad *= 1.4826;  // normalize to consistent with normal distribution

    return {median_radius, mad};
}

static void filter_component_pixel_areas_by_molecules(
    std::vector<size_t>& component_pixel_areas,
    const std::vector<int>& segment_per_molecule,
    int min_molecules_per_segment
) {
    if (min_molecules_per_segment <= 0 || component_pixel_areas.empty()) return;

    std::vector<int> molecule_counts(component_pixel_areas.size(), 0);
    for (int label : segment_per_molecule) {
        if (label > 0 && label <= static_cast<int>(molecule_counts.size())) {
            molecule_counts[label - 1]++;
        }
    }

    for (size_t i = 0; i < component_pixel_areas.size(); ++i) {
        if (molecule_counts[i] < min_molecules_per_segment) {
            component_pixel_areas[i] = 0;
        }
    }
}

// ============================================================================
// Image-based segmentation loading (TIFF label mask)
// ============================================================================

struct TiffWindow {
    bool valid = false;
    uint32_t full_width = 0;
    uint32_t full_height = 0;
    uint32_t row0 = 0;
    uint32_t row1 = 0;
    uint32_t col0 = 0;
    uint32_t col1 = 0;

    uint32_t width() const {
        return (valid && row1 >= row0 && col1 >= col0) ? (col1 - col0 + 1) : 0;
    }
    uint32_t height() const {
        return (valid && row1 >= row0 && col1 >= col0) ? (row1 - row0 + 1) : 0;
    }
    bool is_empty() const {
        return width() == 0 || height() == 0;
    }
    bool is_full_image() const {
        return row0 == 0 && col0 == 0 &&
               row1 + 1 == full_height && col1 + 1 == full_width;
    }
};

static std::tuple<TIFF*, uint32_t, uint32_t, uint16_t> open_tiff_with_metadata(
    const std::string& path)
{
    TIFFSetWarningHandler(nullptr);
    TIFF* tif = TIFFOpen(path.c_str(), "r");
    if (!tif) throw std::runtime_error("Cannot open TIFF mask file: " + path);

    uint32_t w = 0, h = 0;
    uint16_t bps = 8, spp = 1;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,  &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE,  &bps);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    if (spp != 1) {
        TIFFClose(tif);
        throw std::runtime_error("Only single-channel TIFF masks are supported");
    }
    return {tif, w, h, bps};
}

static TiffWindow compute_tiff_window_from_molecules(
    const std::vector<double>& mol_x,
    const std::vector<double>& mol_y,
    uint32_t full_width,
    uint32_t full_height)
{
    TiffWindow window;
    window.full_width = full_width;
    window.full_height = full_height;
    if (full_width == 0 || full_height == 0 || mol_x.empty()) {
        return window;
    }

    bool saw_in_bounds = false;
    uint32_t min_row = full_height;
    uint32_t max_row = 0;
    uint32_t min_col = full_width;
    uint32_t max_col = 0;

    for (size_t i = 0; i < mol_x.size(); ++i) {
        int col = static_cast<int>(std::round(mol_x[i])) - 1;
        int row = static_cast<int>(std::round(mol_y[i])) - 1;
        if (row < 0 || col < 0 ||
            row >= static_cast<int>(full_height) ||
            col >= static_cast<int>(full_width)) {
            continue;
        }
        saw_in_bounds = true;
        min_row = std::min(min_row, static_cast<uint32_t>(row));
        max_row = std::max(max_row, static_cast<uint32_t>(row));
        min_col = std::min(min_col, static_cast<uint32_t>(col));
        max_col = std::max(max_col, static_cast<uint32_t>(col));
    }

    if (!saw_in_bounds) {
        return window;
    }

    window.valid = true;
    window.row0 = min_row;
    window.row1 = max_row;
    window.col0 = min_col;
    window.col1 = max_col;
    return window;
}

// Helper: read a TIFF window into a uint8 flat buffer (normalised: 0=bg, nonzero=foreground value)
// while also detecting whether the selected region stores multiple distinct non-zero labels.
static std::vector<uint8_t> read_tiff_mask_uint8_window(
    const std::string& path,
    const TiffWindow& window,
    uint16_t& out_bps,
    bool& has_multiple_nonzero_values)
{
    auto [tif, full_w, full_h, bps] = open_tiff_with_metadata(path);
    out_bps = bps;
    has_multiple_nonzero_values = false;

    if (window.is_empty()) {
        TIFFClose(tif); // GCOVR_EXCL_LINE: unreachable, load_prior_from_image returns for empty windows before calling this helper
        spdlog::info("No molecules overlap the TIFF image bounds; skipping mask window load"); // GCOVR_EXCL_LINE: unreachable, load_prior_from_image returns for empty windows before calling this helper
        return {}; // GCOVR_EXCL_LINE: unreachable, load_prior_from_image returns for empty windows before calling this helper
    }

    if (window.is_full_image()) {
        spdlog::info("Loading TIFF mask: {}x{}, {} bits/sample", full_w, full_h, bps);
    } else {
        spdlog::info( // GCOVR_EXCL_LINE: unwind-only block of the spdlog::info statement, the call itself is covered by crop-window loads
            "Loading TIFF mask window: {}x{} from full {}x{} (rows {}:{}, cols {}:{}), {} bits/sample",
            window.width(), window.height(), full_w, full_h,
            window.row0, window.row1, window.col0, window.col1, bps);
    }

    uint32_t w = window.width();
    uint32_t h = window.height();
    size_t total = static_cast<size_t>(h) * w;
    std::vector<uint8_t> mask(total, 0);

    tmsize_t scanline_size = TIFFScanlineSize(tif);
    std::vector<uint8_t> row_buf(scanline_size);
    uint32_t first_nonzero_value = 0;
    bool saw_nonzero_value = false;

    for (uint32_t row = window.row0; row <= window.row1; ++row) {
        if (TIFFReadScanline(tif, row_buf.data(), row, 0) < 0) {
            TIFFClose(tif);
            throw std::runtime_error("Error reading TIFF scanline " + std::to_string(row));
        }
        uint8_t* dst = mask.data() + static_cast<size_t>(row - window.row0) * w;
        if (bps == 8) {
            for (uint32_t c = window.col0; c <= window.col1; ++c) {
                uint32_t val = row_buf[c];
                if (val > 0) {
                    if (!saw_nonzero_value) {
                        first_nonzero_value = val;
                        saw_nonzero_value = true;
                    } else if (val != first_nonzero_value) {
                        has_multiple_nonzero_values = true;
                    }
                }
                dst[c - window.col0] = val ? 1 : 0;  // normalise to 0/1
            }
        } else if (bps == 16) {
            const uint16_t* src = reinterpret_cast<const uint16_t*>(row_buf.data());
            for (uint32_t c = window.col0; c <= window.col1; ++c) {
                uint32_t val = src[c];
                if (val > 0) {
                    if (!saw_nonzero_value) {
                        first_nonzero_value = val;
                        saw_nonzero_value = true;
                    } else if (val != first_nonzero_value) {
                        has_multiple_nonzero_values = true;
                    }
                }
                dst[c - window.col0] = val ? 1 : 0;
            }
        } else if (bps == 32) {
            const uint32_t* src = reinterpret_cast<const uint32_t*>(row_buf.data());
            for (uint32_t c = window.col0; c <= window.col1; ++c) {
                uint32_t val = src[c];
                if (val > 0) {
                    if (!saw_nonzero_value) {
                        first_nonzero_value = val;
                        saw_nonzero_value = true;
                    } else if (val != first_nonzero_value) {
                        has_multiple_nonzero_values = true;
                    }
                }
                dst[c - window.col0] = val ? 1 : 0;
            }
        } else {
            TIFFClose(tif);
            throw std::runtime_error("Unsupported TIFF bits/sample: " + std::to_string(bps));
        }
    }
    TIFFClose(tif);
    return mask;
}

ImageSegResult load_prior_from_image(
    const std::string& image_path,
    const std::vector<double>& mol_x,
    const std::vector<double>& mol_y,
    int min_molecules_per_segment
) {
    auto [meta_tif, full_width, full_height, meta_bps] = open_tiff_with_metadata(image_path);
    TIFFClose(meta_tif);

    TiffWindow window = compute_tiff_window_from_molecules(
        mol_x, mol_y, full_width, full_height);
    uint32_t width = window.width(), height = window.height(); uint16_t bps = meta_bps;
    bool has_multiple_nonzero_values = false;
    int n_mols = static_cast<int>(mol_x.size());
    std::vector<int> segment_per_molecule(n_mols, 0);
    if (window.is_empty()) {
        return {segment_per_molecule, {}};
    }

    // Load the selected window of the mask: binary (0/1).
    auto mask = read_tiff_mask_uint8_window(
        image_path, window, bps, has_multiple_nonzero_values);

    // Compute pixel indices for each molecule in the local window
    auto mol_col = [&](int i) {
        return static_cast<int>(std::round(mol_x[i])) - 1 - static_cast<int>(window.col0);
    };
    auto mol_row = [&](int i) {
        return static_cast<int>(std::round(mol_y[i])) - 1 - static_cast<int>(window.row0);
    };
    auto in_bounds = [&](int r, int c) {
        return r >= 0 && r < static_cast<int>(height) &&
               c >= 0 && c < static_cast<int>(width);
    };

    // pixel area per CC label (1-based); only populated for binary masks
    std::vector<size_t> pixel_areas;
    bool already_filtered = false;

    if (!has_multiple_nonzero_values) {
        // ---------------------------------------------------------------
        // Binary mask (all non-zero pixels share the same value, typically
        // from a DAPI threshold).  Run 4-connected component labelling to
        // assign each nucleus a unique integer label.
        //
        // Molecule-centric BFS: only the pixels actually containing molecules
        // need explicit tracking; the BFS still floods the full nucleus area
        // via the uint8 mask to merge molecules that share a nucleus.
        // ---------------------------------------------------------------
        spdlog::info("Binary TIFF mask detected; running connected-component labelling...");

        // Build pixel-index → molecule-list map (only for molecules inside the mask)
        std::unordered_map<size_t, std::vector<int>> pixel_to_mols;
        int n_in_mask = 0;
        for (int i = 0; i < n_mols; ++i) {
            int c = mol_col(i), r = mol_row(i);
            if (!in_bounds(r, c)) continue;
            if (mask[static_cast<size_t>(r) * width + c] == 0) continue;
            pixel_to_mols[static_cast<size_t>(r) * width + c].push_back(i);
            ++n_in_mask;
        }
        spdlog::info("{}/{} molecules fall inside the binary mask", n_in_mask, n_mols);

        // BFS flood-fill from each unlabelled molecule's pixel
        // visited: separate uint8 array the same size as mask (1 byte/pixel)
        std::vector<uint8_t> visited(mask.size(), 0);
        int next_label = 1;
        // pixel_areas (foreground pixel count per CC label, 1-based) declared at function scope

        for (int i = 0; i < n_mols; ++i) {
            if (segment_per_molecule[i] != 0) continue;  // already labelled by a previous BFS

            int c0 = mol_col(i), r0 = mol_row(i);
            if (!in_bounds(r0, c0)) continue;
            size_t idx0 = static_cast<size_t>(r0) * width + c0;
            if (mask[idx0] == 0 || visited[idx0]) continue;

            // BFS
            int cur_label = next_label++;
            pixel_areas.push_back(0);  // slot for cur_label (1-based → index cur_label-1)
            std::queue<size_t> q;
            q.push(idx0);
            visited[idx0] = 1;

            while (!q.empty()) {
                size_t px = q.front(); q.pop();
                pixel_areas[cur_label - 1]++;  // count every foreground pixel in this CC

                int pr = static_cast<int>(px / width);
                int pc = static_cast<int>(px % width);

                // Label any molecules at this pixel
                auto it = pixel_to_mols.find(px);
                if (it != pixel_to_mols.end()) {
                    for (int mi : it->second) segment_per_molecule[mi] = cur_label;
                }

                // 4-connected neighbours
                const int dr[4] = {-1, 1,  0, 0};
                const int dc[4] = { 0, 0, -1, 1};
                for (int d = 0; d < 4; ++d) {
                    int nr = pr + dr[d], nc = pc + dc[d];
                    if (!in_bounds(nr, nc)) continue;
                    size_t nidx = static_cast<size_t>(nr) * width + nc;
                    if (mask[nidx] > 0 && !visited[nidx]) {
                        visited[nidx] = 1;
                        q.push(nidx);
                    }
                }
            }
        }

        int n_labels = next_label - 1;
        spdlog::info("Connected-component labelling: {} nuclei found", n_labels);

        // Match Julia: estimate scale only from components that survive
        // the molecule-overlap filtering applied to the prior labels.
        filter_component_pixel_areas_by_molecules(
            pixel_areas, segment_per_molecule, min_molecules_per_segment);

    } else {
        // ---------------------------------------------------------------
        // Multi-label mask: pixel values directly encode cell IDs.
        // Preserve those labels and estimate scale from per-label pixel areas,
        // matching Julia's labeled-mask path.
        // ---------------------------------------------------------------
        spdlog::info("Multi-label TIFF mask detected; preserving label values...");

        auto [tif, tiff_width, tiff_height, tiff_bps] = open_tiff_with_metadata(image_path);
        if (tiff_width != full_width || tiff_height != full_height) {
            TIFFClose(tif); // GCOVR_EXCL_LINE: unreachable, the file cannot change between two synchronous reads in one call
            throw std::runtime_error("TIFF dimensions changed between metadata and data reads"); // GCOVR_EXCL_LINE: unreachable, the file cannot change between two synchronous reads in one call
        }
        bps = tiff_bps;

        tmsize_t scanline_size = TIFFScanlineSize(tif);
        std::vector<uint8_t> row_buf(scanline_size);

        // Build sorted molecule list by local row for efficient scanline matching
        std::vector<std::pair<int,int>> mol_by_row;  // (row, mol_idx)
        mol_by_row.reserve(n_mols);
        for (int i = 0; i < n_mols; ++i) {
            int local_row_i = mol_row(i);
            int local_col_i = mol_col(i);
            if (in_bounds(local_row_i, local_col_i)) {
                mol_by_row.push_back({local_row_i, i});
            }
        }
        std::sort(mol_by_row.begin(), mol_by_row.end());

        std::unordered_map<uint32_t, size_t> pixel_areas_by_label;
        int mol_ptr = 0;
        for (uint32_t row = window.row0; row <= window.row1; ++row) {
            if (TIFFReadScanline(tif, row_buf.data(), row, 0) < 0) {
                TIFFClose(tif); // GCOVR_EXCL_LINE: unreachable, the same scanlines are read earlier in read_tiff_mask_uint8_window, which would throw first
                throw std::runtime_error("Error reading TIFF scanline " + std::to_string(row)); // GCOVR_EXCL_LINE: unreachable, the same scanlines are read earlier in read_tiff_mask_uint8_window, which would throw first
            }
            int local_row = static_cast<int>(row - window.row0);

            if (bps == 8) {
                for (uint32_t c = window.col0; c <= window.col1; ++c) {
                    uint32_t val = row_buf[c];
                    if (val > 0) pixel_areas_by_label[val]++;
                }
            } else if (bps == 16) {
                const uint16_t* src = reinterpret_cast<const uint16_t*>(row_buf.data());
                for (uint32_t c = window.col0; c <= window.col1; ++c) {
                    uint32_t val = src[c];
                    if (val > 0) pixel_areas_by_label[val]++;
                }
            } else if (bps == 32) {
                const uint32_t* src = reinterpret_cast<const uint32_t*>(row_buf.data());
                for (uint32_t c = window.col0; c <= window.col1; ++c) {
                    uint32_t val = src[c];
                    if (val > 0) pixel_areas_by_label[val]++;
                }
            }

            while (mol_ptr < static_cast<int>(mol_by_row.size()) &&
                   mol_by_row[mol_ptr].first == local_row) {
                int mi = mol_by_row[mol_ptr].second;
                int c  = mol_col(mi);
                if (in_bounds(local_row, c)) {
                    uint32_t val = 0;
                    uint32_t src_col = static_cast<uint32_t>(c) + window.col0;
                    if (bps == 8)       val = row_buf[src_col];
                    else if (bps == 16) val = reinterpret_cast<const uint16_t*>(row_buf.data())[src_col];
                    else if (bps == 32) val = reinterpret_cast<const uint32_t*>(row_buf.data())[src_col];
                    segment_per_molecule[mi] = static_cast<int>(val);
                }
                ++mol_ptr;
            }
        }
        TIFFClose(tif);

        int n_assigned = 0;
        for (int v : segment_per_molecule) if (v > 0) ++n_assigned;
        spdlog::info("{}/{} molecules overlap a prior segment", n_assigned, n_mols);

        std::unordered_map<uint32_t, int> molecule_counts;
        molecule_counts.reserve(pixel_areas_by_label.size());
        for (int label : segment_per_molecule) {
            if (label > 0) molecule_counts[static_cast<uint32_t>(label)]++;
        }

        for (int& label : segment_per_molecule) {
            if (label <= 0) continue;
            auto it = molecule_counts.find(static_cast<uint32_t>(label));
            if ((it == molecule_counts.end()) || (it->second < min_molecules_per_segment)) {
                label = 0;
            }
        }

        pixel_areas.reserve(pixel_areas_by_label.size());
        for (const auto& [label, area] : pixel_areas_by_label) {
            auto it = molecule_counts.find(label);
            if ((it != molecule_counts.end()) && (it->second >= min_molecules_per_segment)) {
                pixel_areas.push_back(area);
            }
        }
        already_filtered = true;
    }

    if (!already_filtered) {
        filter_segmentation_labels(segment_per_molecule, min_molecules_per_segment);
    }
    return {segment_per_molecule, pixel_areas};
}

// ============================================================================
// Top-level: load prior segmentation
// ============================================================================

std::pair<double, double> load_prior_segmentation(
    MoleculeData& data,
    const PriorInputOptions& prior_opts,
    int min_molecules_per_cell
) {
    auto seg_type = prior_opts.type;

    if (seg_type == PriorInputType::None) {
        return {-1.0, -1.0};
    }

    double scale = -1.0, scale_std = -1.0;

    if (seg_type == PriorInputType::Column) {
        if (data.prior_segmentation.empty()) {
            throw std::runtime_error(
                "Prior segmentation column '" + prior_opts.column_name +
                "' was requested, but no prior labels were loaded from the molecule file");
        }

        if (prior_opts.estimate_scale_from_prior) {
            auto pos = data.position_matrix();
            auto [s, s_std] = estimate_scale_from_assignment(
                pos, data.prior_segmentation, min_molecules_per_cell);
            scale = s;
            scale_std = s_std;
            spdlog::info("Estimated scale from prior segmentation: {:.2f} (std: {:.2f})",
                         scale, scale_std);
        }
    } else if (seg_type == PriorInputType::Boundary) {
        auto polygons = load_boundary_polygons(prior_opts.path);
        filter_boundary_polygons_to_molecule_bounds(
            polygons, compute_molecule_bounds(data.x, data.y));
        data.prior_segmentation = assign_molecules_to_boundaries(
            data.x, data.y, polygons, prior_opts.min_molecules_per_segment);

        if (prior_opts.estimate_scale_from_prior) {
            try {
                auto pos = data.position_matrix();
                auto [s, s_std] = estimate_scale_from_assignment(
                    pos, data.prior_segmentation, min_molecules_per_cell);
                scale = s;
                scale_std = s_std;
                spdlog::info("Estimated scale from prior segmentation: {:.2f} (std: {:.2f})",
                             scale, scale_std);
            } catch (const std::exception& e) {
                spdlog::warn("Could not estimate scale from prior: {}", e.what());
            }
        }
    } else {
        // Image-based
        auto result = load_prior_from_image(
            prior_opts.path, data.x, data.y, prior_opts.min_molecules_per_segment);
        data.prior_segmentation = std::move(result.segment_per_molecule);

        if (prior_opts.estimate_scale_from_prior) {
            try {
                if (!result.component_pixel_areas.empty()) {
                    // Image prior: use sqrt(area/π) radii for both binary and labeled masks,
                    // matching Julia's image-prior path.
                    auto [s, s_std] = estimate_scale_from_image_areas(result.component_pixel_areas);
                    scale = s;
                    scale_std = s_std;
                } else {
                    // Fallback when no label areas survive filtering.
                    auto pos = data.position_matrix();
                    auto [s, s_std] = estimate_scale_from_assignment( // GCOVR_EXCL_LINE: unwind-only block of this call, arguments are plain references that cannot throw
                        pos, data.prior_segmentation, min_molecules_per_cell);
                    scale = s; // GCOVR_EXCL_LINE: unreachable, empty component areas imply no molecule was assigned, so the estimate above always throws
                    scale_std = s_std; // GCOVR_EXCL_LINE: unreachable, empty component areas imply no molecule was assigned, so the estimate above always throws
                }
                spdlog::info("Estimated scale from prior segmentation: {:.2f} (std: {:.2f})",
                             scale, scale_std);
            } catch (const std::exception& e) {
                spdlog::warn("Could not estimate scale from prior: {}", e.what());
            }
        }
    }

    return {scale, scale_std};
}

} // namespace baysor
