#include "baysor/processing/data_processing/boundary_estimation.h"
#include "baysor/processing/data_processing/boundary_estimation_internal.h"

#include "baysor/processing/data_processing/triangulation.h"
#include "baysor/processing/utils/utils.h"
#include "baysor/utils/general.h"
#include "baysor/utils/thread_pool.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <spdlog/spdlog.h>
#include <sstream>
#include <unordered_map>

namespace baysor {

namespace {

using BBox = std::array<double, 4>; // xmin, xmax, ymin, ymax
using Edge = std::pair<int, int>;
using Triangle = std::array<int, 3>;

using CgalKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using CgalVertexBase = CGAL::Triangulation_vertex_base_with_info_2<int, CgalKernel>;
using CgalFaceBase = CGAL::Triangulation_face_base_2<CgalKernel>;
using CgalTds = CGAL::Triangulation_data_structure_2<CgalVertexBase, CgalFaceBase>;
using CgalDelaunay = CGAL::Delaunay_triangulation_2<CgalKernel, CgalTds>;
using CgalPoint = CgalKernel::Point_2;
using CgalPointWithInfo = std::pair<CgalPoint, int>;

inline std::uint64_t edge_key(int a, int b) {
    const int lo = std::min(a, b);
    const int hi = std::max(a, b);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32)
         | static_cast<std::uint32_t>(hi);
}

inline Edge canonical_edge(int a, int b) {
    return {std::min(a, b), std::max(a, b)};
}

std::string default_cell_name(int cid) {
    return std::to_string(cid);
}

std::string format_double(double x) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(6);
    oss << x;
    std::string s = oss.str();
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    return s.empty() ? "0" : s;
}

double quantile_vec(std::vector<double> vals, double p) {
    if (vals.empty()) return 0.0;
    std::sort(vals.begin(), vals.end());
    const double idx = p * static_cast<double>(vals.size() - 1);
    const int lo = static_cast<int>(std::floor(idx));
    const int hi = static_cast<int>(std::ceil(idx));
    if (lo == hi) return vals[lo];
    const double frac = idx - lo;
    return vals[lo] * (1.0 - frac) + vals[hi] * frac;
}

Eigen::MatrixXd subset_columns(const Eigen::MatrixXd& mat, const std::vector<int>& ids) {
    Eigen::MatrixXd out(mat.rows(), static_cast<int>(ids.size()));
    for (int i = 0; i < static_cast<int>(ids.size()); ++i) {
        out.col(i) = mat.col(ids[i]);
    }
    return out;
}

std::vector<BBox> get_boundary_box_per_cell(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& assignment,
    int max_label
) {
    auto ids_per_cell = split_ids(assignment, max_label, /*drop_zero=*/true);
    std::vector<BBox> bboxes(max_label, {
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()
    });

    for (int cid = 0; cid < max_label; ++cid) {
        const auto& ids = ids_per_cell[cid];
        if (ids.empty()) continue;
        for (int idx : ids) {
            const double x = pos_data(0, idx);
            const double y = pos_data(1, idx);
            bboxes[cid][0] = std::min(bboxes[cid][0], x);
            bboxes[cid][1] = std::max(bboxes[cid][1], x);
            bboxes[cid][2] = std::min(bboxes[cid][2], y);
            bboxes[cid][3] = std::max(bboxes[cid][3], y);
        }
    }

    return bboxes;
}

std::vector<std::vector<int>> extract_ids_per_bbox(
    const Eigen::VectorXd& xvals,
    const Eigen::VectorXd& yvals,
    const std::vector<BBox>& bboxes
) {
    const int n = static_cast<int>(xvals.size());
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return xvals[a] < xvals[b];
    });

    std::vector<std::vector<int>> ids_per_cell(bboxes.size());
    for (int cid = 0; cid < static_cast<int>(bboxes.size()); ++cid) {
        const auto& bb = bboxes[cid];
        if (!std::isfinite(bb[0])) continue;

        const auto it = std::lower_bound(order.begin(), order.end(), bb[0], [&](int idx, double x) {
            return xvals[idx] < x;
        });
        for (auto jt = it; jt != order.end(); ++jt) {
            const int idx = *jt;
            const double x = xvals[idx];
            if (x > bb[1]) break;
            const double y = yvals[idx];
            if (y >= bb[2] && y <= bb[3]) {
                ids_per_cell[cid].push_back(idx);
            }
        }
    }

    return ids_per_cell;
}

int get_n_points_in_triangle(const Eigen::MatrixXd& points, const Eigen::Matrix<double, 2, 3>& tri) {
    if (points.cols() == 0) return 0;

    const Eigen::Vector2d a = tri.col(0);
    const Eigen::Vector2d v0 = tri.col(2) - a;
    const Eigen::Vector2d v1 = tri.col(1) - a;

    const double dot00 = v0.dot(v0);
    const double dot01 = v0.dot(v1);
    const double dot11 = v1.dot(v1);
    const double denom = dot00 * dot11 - dot01 * dot01;
    if (std::abs(denom) < 1e-12) return 0;
    const double inv_denom = 1.0 / denom;

    int n_inner = 0;
    for (int i = 0; i < points.cols(); ++i) {
        const Eigen::Vector2d v2 = points.col(i) - a;
        const double dot02 = v0.dot(v2);
        const double dot12 = v1.dot(v2);
        const double u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
        const double v = (dot00 * dot12 - dot01 * dot02) * inv_denom;
        if ((u >= 0.0) && (v >= 0.0) && (u + v < 1.0)) {
            ++n_inner;
        }
    }
    return n_inner;
}

std::vector<Triangle> extract_triangle_verts(CgalDelaunay& dt) {
    std::vector<Triangle> triangles;
    triangles.reserve(static_cast<size_t>(dt.number_of_faces()));
    for (auto fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) {
        triangles.push_back({ // GCOVR_EXCL_LINE: dead GCC block; push_back call counted on lines 179-181
            fit->vertex(0)->info(),
            fit->vertex(1)->info(),
            fit->vertex(2)->info()
        });
    }
    return triangles;
}

std::vector<Edge> extract_border_edges(const std::vector<Triangle>& triangles) {
    std::unordered_map<std::uint64_t, int> edge_counts;
    edge_counts.reserve(triangles.size() * 6 + 1);

    std::vector<Edge> tri_edges;
    tri_edges.reserve(triangles.size() * 3);
    for (const auto& tri : triangles) {
        const Edge e1 = canonical_edge(tri[0], tri[1]);
        const Edge e2 = canonical_edge(tri[1], tri[2]);
        const Edge e3 = canonical_edge(tri[2], tri[0]);
        tri_edges.push_back(e1);
        tri_edges.push_back(e2);
        tri_edges.push_back(e3);
        edge_counts[edge_key(e1.first, e1.second)]++;
        edge_counts[edge_key(e2.first, e2.second)]++;
        edge_counts[edge_key(e3.first, e3.second)]++;
    }

    std::vector<Edge> border_edges;
    border_edges.reserve(tri_edges.size());
    for (const auto& e : tri_edges) {
        if (edge_counts[edge_key(e.first, e.second)] == 1) {
            border_edges.push_back(e);
        }
    }
    return border_edges;
}

std::vector<Edge> find_border_without_admixture(
    const std::vector<Triangle>& triangles,
    const Eigen::MatrixXd& pos_data,
    const Eigen::MatrixXd& non_cell_pos,
    int max_iters
) {
    std::vector<std::array<Edge, 3>> edges_per_tri(triangles.size());
    std::unordered_map<std::uint64_t, int> edge_counts;
    edge_counts.reserve(triangles.size() * 6 + 1);

    for (int i = 0; i < static_cast<int>(triangles.size()); ++i) {
        const auto& tri = triangles[i];
        edges_per_tri[i] = {
            canonical_edge(tri[0], tri[1]),
            canonical_edge(tri[1], tri[2]),
            canonical_edge(tri[2], tri[0])
        };
        for (const auto& e : edges_per_tri[i]) {
            edge_counts[edge_key(e.first, e.second)]++;
        }
    }

    std::unordered_map<int, int> n_borders_per_node;
    for (const auto& [key, n] : edge_counts) {
        if (n != 1) continue;
        const int a = static_cast<int>(key >> 32);
        const int b = static_cast<int>(key & 0xffffffffu);
        n_borders_per_node[a] = n_borders_per_node[a] + 1;
        n_borders_per_node[b] = n_borders_per_node[b] + 1;
    }

    std::vector<char> excluded(triangles.size(), 0);
    for (int iter = 0; iter < max_iters; ++iter) {
        bool converged = true;
        for (int i = 0; i < static_cast<int>(triangles.size()); ++i) {
            if (excluded[i]) continue;

            int n_borders = 0;
            for (const auto& e : edges_per_tri[i]) {
                n_borders += (edge_counts[edge_key(e.first, e.second)] == 1);
            }
            if (n_borders == 0) continue;
            if (n_borders > 1) continue;

            Eigen::Matrix<double, 2, 3> tri_pts;
            tri_pts.col(0) = pos_data.col(triangles[i][0]);
            tri_pts.col(1) = pos_data.col(triangles[i][1]);
            tri_pts.col(2) = pos_data.col(triangles[i][2]);
            const int n_admix = get_n_points_in_triangle(non_cell_pos, tri_pts);
            if (n_admix <= 0) continue;

            bool skip = false;
            for (const auto& e : edges_per_tri[i]) {
                const auto key = edge_key(e.first, e.second);
                if (edge_counts[key] != 2) continue;
                const bool both_have_other_border_edges =
                    (n_borders_per_node[e.first] == 2) && (n_borders_per_node[e.second] == 2);
                if (both_have_other_border_edges) {
                    skip = true;
                    break;
                }
            }
            if (skip) continue;

            converged = false;
            excluded[i] = 1;
            for (const auto& e : edges_per_tri[i]) {
                const auto key = edge_key(e.first, e.second);
                edge_counts[key] -= 1;
                const int new_count = edge_counts[key];
                if (new_count == 0) {
                    n_borders_per_node[e.first] -= 1;
                    n_borders_per_node[e.second] -= 1;
                } else if (new_count == 1) {
                    n_borders_per_node[e.first] += 1;
                    n_borders_per_node[e.second] += 1;
                }
            }
        }

        if (converged) break;
        if (iter == max_iters - 1) {
            spdlog::warn("Polygon filtering did not converge within {} iterations", max_iters);
        }
    }

    std::vector<Edge> border_edges;
    border_edges.reserve(edge_counts.size());
    for (const auto& [key, n] : edge_counts) {
        if (n != 1) continue;
        border_edges.push_back({
            static_cast<int>(key >> 32), // GCOVR_EXCL_LINE: dead GCC block; expression counted on line 307
            static_cast<int>(key & 0xffffffffu)
        });
    }
    return border_edges;
}

std::vector<int> border_edges_to_poly(const std::vector<Edge>& border_edges, int max_border_len) {
    if (border_edges.size() <= 2) return {};

    std::unordered_map<int, std::pair<int, int>> adjacency;
    adjacency.reserve(border_edges.size() * 2 + 1);

    for (const auto& [a, b] : border_edges) {
        auto ita = adjacency.find(a);
        if (ita == adjacency.end()) {
            adjacency[a] = {b, -1};
        } else {
            ita->second.second = b;
        }

        auto itb = adjacency.find(b);
        if (itb == adjacency.end()) {
            adjacency[b] = {a, -1};
        } else {
            itb->second.second = a;
        }
    }

    if (adjacency.empty()) return {};
    int start = adjacency.begin()->first;
    int cur = start;
    int prev = -1;

    std::vector<int> poly;
    poly.reserve(adjacency.size());
    poly.push_back(start);

    for (int iter = 0; iter < max_border_len; ++iter) {
        const auto it = adjacency.find(cur);
        if (it == adjacency.end()) return {};
        const auto [a, b] = it->second;
        const int next = (a == prev) ? b : a;
        if (next < 0) return {};
        prev = cur;
        cur = next;
        if (cur == start) {
            return poly;
        }
        poly.push_back(cur);
    }

    spdlog::warn("Could not build a polygon border of size {}", adjacency.size());
    return {};
}

PolygonCollection build_polygons_for_cells(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& cell_labels,
    const std::vector<std::string>* cell_names,
    double offset_rel
) {
    const int n = static_cast<int>(pos_data.cols());
    if (n < 3) return {};

    const int max_label = count_array(cell_labels, -1, /*drop_zero=*/true).size();
    if (max_label <= 0) return {};

    const Eigen::MatrixXd pos2d = pos_data.topRows(2);
    const Eigen::MatrixXd norm_pts = normalize_points(pos2d);

    auto bboxes = get_boundary_box_per_cell(pos2d, cell_labels, max_label);
    auto ids_per_bbox = extract_ids_per_bbox(pos2d.row(0).transpose(), pos2d.row(1).transpose(), bboxes);
    auto mids_per_cell = split_ids(cell_labels, max_label, /*drop_zero=*/true);

    auto knn = knn_parallel(pos2d, pos2d, 2, true);
    double mean_nn_dist = 1.0;
    if (knn.n > 0) {
        double sum = 0.0;
        int count = 0;
        for (int i = 0; i < knn.n; ++i) {
            if (knn.k >= 2) {
                sum += knn.dist_row(i)[1];
                ++count;
            }
        }
        if (count > 0) mean_nn_dist = sum / count;
    }
    const double offset = mean_nn_dist * offset_rel;

    std::vector<Eigen::MatrixXd> polygons_by_cell(max_label);

    // Fallback polygon for a cell whose free-form boundary estimation failed
    // or produced fewer than 3 vertices: the axis-aligned bounding box of the
    // cell's molecules, padded by the local nearest-neighbour offset. Every
    // cell with at least one assigned molecule must have a polygon, otherwise
    // the segmentation CSV and the polygons JSON disagree and Xenium Ranger
    // rejects the import (kharchenkolab/Baysor#165).
    auto fallback_polygon = [&](int cid) {
        const auto& cell_ids = mids_per_cell[cid - 1];
        double min_x = std::numeric_limits<double>::infinity();
        double max_x = -std::numeric_limits<double>::infinity();
        double min_y = std::numeric_limits<double>::infinity();
        double max_y = -std::numeric_limits<double>::infinity();
        for (int mol : cell_ids) {
            const double x = pos2d(0, mol);
            const double y = pos2d(1, mol);
            min_x = std::min(min_x, x);
            max_x = std::max(max_x, x);
            min_y = std::min(min_y, y);
            max_y = std::max(max_y, y);
        }
        const double pad = std::max(offset, 1e-6);
        Eigen::MatrixXd poly(2, 4);
        poly.col(0) = Eigen::Vector2d(min_x - pad, min_y - pad);
        poly.col(1) = Eigen::Vector2d(max_x + pad, min_y - pad);
        poly.col(2) = Eigen::Vector2d(max_x + pad, max_y + pad);
        poly.col(3) = Eigen::Vector2d(min_x - pad, max_y + pad);
        return poly;
    };
    auto emit_fallback = [&](int cid) {
        polygons_by_cell[cid - 1] = fallback_polygon(cid);
    };

    // Nested calls inside a parallel region run serially (this replaces the
    // former `omp_in_parallel` guard).
    parallel_for(1, max_label + 1, 32, [&](int cid) {
        const auto& cell_ids = mids_per_cell[cid - 1];
        if (cell_ids.empty()) return;

        Eigen::MatrixXd poly;
        if (cell_ids.size() == 1) {
            const Eigen::Vector2d p = pos2d.col(cell_ids[0]);
            poly.resize(2, 4);
            poly.col(0) = p + Eigen::Vector2d(offset, 0.0);
            poly.col(1) = p + Eigen::Vector2d(-offset, 0.0);
            poly.col(2) = p + Eigen::Vector2d(0.0, offset);
            poly.col(3) = p + Eigen::Vector2d(0.0, -offset);
        } else if (cell_ids.size() == 2) {
            const Eigen::Vector2d p1 = pos2d.col(cell_ids[0]);
            const Eigen::Vector2d p2 = pos2d.col(cell_ids[1]);
            const Eigen::Vector2d center = (p1 + p2) / 2.0;
            poly.resize(2, 4);
            poly.col(0) = p1;
            poly.col(1) = p2;
            poly.col(2) = center + Eigen::Vector2d(offset, 0.0);
            poly.col(3) = center + Eigen::Vector2d(-offset, 0.0);
        } else {
            const auto& mids = ids_per_bbox[cid - 1];
            if (mids.empty()) { emit_fallback(cid); return; }

            const Eigen::MatrixXd bbox_pos = subset_columns(pos2d, mids);
            const Eigen::MatrixXd bbox_norm = subset_columns(norm_pts, mids);

            std::vector<int> bbox_labels(mids.size());
            int n_cell_pts = 0;
            for (int i = 0; i < static_cast<int>(mids.size()); ++i) {
                bbox_labels[i] = cell_labels[mids[i]];
                n_cell_pts += (bbox_labels[i] == cid);
            }
            if (n_cell_pts < 3) { emit_fallback(cid); return; }

            std::vector<CgalPointWithInfo> cell_points;
            cell_points.reserve(n_cell_pts);
            std::vector<int> non_cell_ids;
            non_cell_ids.reserve(mids.size() - n_cell_pts);

            for (int i = 0; i < static_cast<int>(mids.size()); ++i) {
                if (bbox_labels[i] == cid) {
                    cell_points.push_back({CgalPoint(bbox_norm(0, i), bbox_norm(1, i)), i});
                } else {
                    non_cell_ids.push_back(i);
                }
            }

            CgalDelaunay dt;
            dt.insert(cell_points.begin(), cell_points.end());
            auto triangles = extract_triangle_verts(dt);
            if (triangles.empty()) { emit_fallback(cid); return; }

            Eigen::MatrixXd non_cell_pos(2, static_cast<int>(non_cell_ids.size()));
            for (int i = 0; i < static_cast<int>(non_cell_ids.size()); ++i) {
                non_cell_pos.col(i) = bbox_pos.col(non_cell_ids[i]);
            }

            auto border_edges = find_border_without_admixture(triangles, bbox_pos, non_cell_pos,
                                                              internal::kBorderFilterMaxIters);
            auto poly_ids = border_edges_to_poly(border_edges, internal::kMaxBorderLength);
            if (poly_ids.empty()) { emit_fallback(cid); return; }

            poly.resize(2, static_cast<int>(poly_ids.size()));
            for (int i = 0; i < static_cast<int>(poly_ids.size()); ++i) {
                poly.col(i) = bbox_pos.col(poly_ids[i]);
            }
        }

        if (poly.cols() < 3) { emit_fallback(cid); return; }
        polygons_by_cell[cid - 1] = std::move(poly);
    });

    PolygonCollection polygons;
    polygons.reserve(max_label);
    for (int cid = 1; cid <= max_label; ++cid) {
        if (polygons_by_cell[cid - 1].cols() == 0) continue;
        // When the caller supplies names, an out-of-range label still has to
        // map to the same name the molecule CSV uses (`cell_<n>`), never to
        // the bare-integer default, so the two files always agree.
        std::string cell_name;
        if (cell_names) {
            cell_name = (cid - 1 < static_cast<int>(cell_names->size()))
                ? (*cell_names)[cid - 1]
                : ("cell_" + std::to_string(cid));
        } else {
            cell_name = default_cell_name(cid);
        }
        polygons[cell_name] = std::move(polygons_by_cell[cid - 1]);
    }

    return polygons;
}

std::vector<std::vector<Eigen::Vector2d>> grid_borders_per_label(
    const Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic>& grid_labels
) {
    const int n_labels = static_cast<int>(grid_labels.maxCoeff());
    std::vector<std::vector<Eigen::Vector2d>> borders(n_labels);
    const int n_rows = static_cast<int>(grid_labels.rows());
    const int n_cols = static_cast<int>(grid_labels.cols());

    static const int dr[4] = {-1, 0, 1, 0};
    static const int dc[4] = {0, 1, 0, -1};

    for (int row = 0; row < n_rows; ++row) {
        for (int col = 0; col < n_cols; ++col) {
            const uint32_t label = grid_labels(row, col);
            if (label == 0) continue;

            for (int k = 0; k < 4; ++k) {
                const int nr = row + dr[k];
                const int nc = col + dc[k];
                if (nr < 0 || nc < 0 || nr >= n_rows || nc >= n_cols ||
                    grid_labels(nr, nc) != label) {
                    borders[label - 1].push_back(Eigen::Vector2d(col + 1.0, row + 1.0));
                    break;
                }
            }
        }
    }

    return borders;
}

} // namespace

PolygonCollection boundary_polygons(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& cell_labels,
    const std::vector<std::string>* cell_names,
    double offset_rel
) {
    if (pos_data.rows() == 3) {
        return build_polygons_for_cells(pos_data, cell_labels, cell_names, offset_rel);
    }
    if (pos_data.rows() != 2) {
        spdlog::error("Only 2D and 3D data are supported for boundary estimation");
        return {};
    }
    return build_polygons_for_cells(pos_data, cell_labels, cell_names, offset_rel);
}

std::vector<Eigen::MatrixXd> boundary_polygons_from_grid(
    const Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic>& grid_labels,
    int grid_step
) {
    const auto borders = grid_borders_per_label(grid_labels);
    std::vector<Eigen::MatrixXd> polys;
    polys.reserve(borders.size());

    for (const auto& border_pts : borders) {
        if (border_pts.size() < 2) {
            polys.emplace_back();
            continue;
        }

        Eigen::MatrixXd coords(2, static_cast<int>(border_pts.size()));
        std::vector<CgalPointWithInfo> pts;
        pts.reserve(border_pts.size());
        for (int i = 0; i < static_cast<int>(border_pts.size()); ++i) {
            coords.col(i) = border_pts[i];
            pts.push_back({CgalPoint(border_pts[i](0), border_pts[i](1)), i});
        }

        CgalDelaunay dt;
        dt.insert(pts.begin(), pts.end());
        auto triangles = extract_triangle_verts(dt);
        auto border_edges = extract_border_edges(triangles);
        auto poly_ids = border_edges_to_poly(border_edges, internal::kMaxBorderLength);
        if (poly_ids.empty()) {
            polys.emplace_back();
            continue;
        }

        Eigen::MatrixXd poly(2, static_cast<int>(poly_ids.size()));
        for (int i = 0; i < static_cast<int>(poly_ids.size()); ++i) {
            poly.col(i) = coords.col(poly_ids[i]) * static_cast<double>(grid_step);
        }
        polys.push_back(std::move(poly));
    }

    return polys;
}

std::pair<PolygonCollection, PolygonStack> boundary_polygons_auto(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& assignment,
    bool estimate_per_z,
    const std::vector<std::string>* cell_names,
    bool verbose,
    int max_z_slices
) {
    if (verbose) {
        spdlog::info("Estimating boundary polygons...");
    }

    PolygonCollection poly_joined = boundary_polygons(pos_data, assignment, cell_names);
    PolygonStack poly_stack;
    poly_stack.push_back({"2d", poly_joined});

    if (!estimate_per_z || pos_data.rows() == 2) {
        return {poly_joined, poly_stack};
    }

    std::vector<double> z_vals(pos_data.cols());
    for (int i = 0; i < pos_data.cols(); ++i) z_vals[i] = pos_data(2, i);

    std::vector<double> unique_z = z_vals;
    std::sort(unique_z.begin(), unique_z.end());
    unique_z.erase(std::unique(unique_z.begin(), unique_z.end()), unique_z.end());

    std::vector<int> z_bins(z_vals.size(), 0);
    std::vector<std::string> layer_names;

    if (static_cast<int>(unique_z.size()) > max_z_slices) {
        if (verbose) {
            spdlog::warn("Too many z values ({}). Binning z-stack into {} layers for polygon estimation. " // GCOVR_EXCL_LINE: dead GCC block; the call executes and is counted on the guard line above
                         "Increase the [plotting] max_z_slices option to use more layers.",
                         unique_z.size(), max_z_slices);
        }
        const double clip = std::min(1.0 / max_z_slices / 4.0, 0.025);
        const double lo = quantile_vec(z_vals, clip);
        const double hi = quantile_vec(z_vals, 1.0 - clip);

        std::vector<double> breaks(max_z_slices + 1, lo);
        const double step = (hi - lo) / max_z_slices;
        for (int i = 0; i <= max_z_slices; ++i) {
            breaks[i] = lo + step * i;
        }

        layer_names.resize(max_z_slices);
        for (int i = 0; i < max_z_slices; ++i) {
            layer_names[i] = "[" + format_double(breaks[i]) + "," + format_double(breaks[i + 1]) + "]";
        }

        for (int i = 0; i < static_cast<int>(z_vals.size()); ++i) {
            int bin = static_cast<int>(std::upper_bound(breaks.begin() + 1, breaks.end() - 1, z_vals[i]) - (breaks.begin() + 1));
            bin = std::max(0, std::min(max_z_slices - 1, bin));
            z_bins[i] = bin;
        }
    } else {
        layer_names.resize(unique_z.size());
        for (int i = 0; i < static_cast<int>(unique_z.size()); ++i) {
            layer_names[i] = format_double(unique_z[i]);
        }
        for (int i = 0; i < static_cast<int>(z_vals.size()); ++i) {
            z_bins[i] = static_cast<int>(std::lower_bound(unique_z.begin(), unique_z.end(), z_vals[i]) - unique_z.begin());
        }
    }

    for (int layer = 0; layer < static_cast<int>(layer_names.size()); ++layer) {
        std::vector<int> ids;
        ids.reserve(z_bins.size());
        for (int i = 0; i < static_cast<int>(z_bins.size()); ++i) {
            if (z_bins[i] == layer) ids.push_back(i);
        }
        if (ids.empty()) continue;

        Eigen::MatrixXd pos_layer = subset_columns(pos_data, ids);
        std::vector<int> ass_layer(ids.size());
        for (int i = 0; i < static_cast<int>(ids.size()); ++i) {
            ass_layer[i] = assignment[ids[i]];
        }
        poly_stack.push_back({layer_names[layer], boundary_polygons(pos_layer, ass_layer, cell_names)});
    }

    return {poly_joined, poly_stack};
}

// Test seam: thin wrappers over the anonymous-namespace helpers above. They
// exist so unit tests can exercise the max_iters / max_border_len warning
// guards, which production callers never override. Behavior is unchanged.
namespace internal {

std::vector<std::pair<int, int>> find_border_without_admixture(
    const std::vector<std::array<int, 3>>& triangles,
    const Eigen::MatrixXd& pos_data,
    const Eigen::MatrixXd& non_cell_pos,
    int max_iters
) {
    return baysor::find_border_without_admixture(triangles, pos_data, non_cell_pos, max_iters);
}

std::vector<int> border_edges_to_poly(
    const std::vector<std::pair<int, int>>& border_edges,
    int max_border_len
) {
    return baysor::border_edges_to_poly(border_edges, max_border_len);
}

} // namespace internal

} // namespace baysor
