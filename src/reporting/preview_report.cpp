#include "baysor/reporting/preview_report.h"
#include "baysor/utils/general.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <third_party/stb_image_write.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <random>
#include <sstream>
#include <vector>
#include <omp.h>

namespace baysor {

static constexpr double PI = 3.14159265358979323846;

struct RasterViewport {
    double xmin, xmax, ymin, ymax, xrange, yrange;
    int width_px, height_px;
};

static inline double normal_pdf(double x, double mu, double sigma) {
    if (sigma <= 0.0) return 0.0;
    double z = (x - mu) / sigma;
    return std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * PI));
}

// ============================================================================
// PNG scatter renderer
// ============================================================================

// Parse a "#RRGGBB" hex string to R,G,B bytes.
static void hex_to_rgb(const std::string& hex, uint8_t& r, uint8_t& g, uint8_t& b) {
    auto h = [&](int pos) {
        char c = hex[pos];
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'A' && c <= 'F') return c - 'A' + 10;
        if (c >= 'a' && c <= 'f') return c - 'a' + 10;
        return 0;
    };
    r = static_cast<uint8_t>((h(1) << 4) | h(2));
    g = static_cast<uint8_t>((h(3) << 4) | h(4));
    b = static_cast<uint8_t>((h(5) << 4) | h(6));
}

// Encode raw RGB pixels to PNG and then to a base64 data URI string.
static std::string pixels_to_base64_png(const std::vector<uint8_t>& pixels,
                                         int width, int height) {
    // Write PNG to memory via stb callback.
    std::vector<uint8_t> png_buf;
    stbi_write_png_compression_level = 6;
    stbi_write_png_to_func(
        [](void* ctx, void* data, int size) { // GCOVR_EXCL_LINE: gcov: lambda entry line carries only an exception-unwind block; the executed entry is reported on the call line
            auto* buf = reinterpret_cast<std::vector<uint8_t>*>(ctx);
            const uint8_t* p = reinterpret_cast<const uint8_t*>(data);
            buf->insert(buf->end(), p, p + size);
        },
        &png_buf, width, height, 3 /*RGB*/, pixels.data(), width * 3);

    // Base64 encode.
    static const char b64[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve(4 * ((png_buf.size() + 2) / 3) + 30);
    out = "data:image/png;base64,";
    size_t i = 0;
    for (; i + 2 < png_buf.size(); i += 3) {
        uint32_t v = (png_buf[i] << 16) | (png_buf[i+1] << 8) | png_buf[i+2];
        out += b64[(v >> 18) & 63];
        out += b64[(v >> 12) & 63];
        out += b64[(v >>  6) & 63];
        out += b64[(v      ) & 63];
    }
    if (i < png_buf.size()) {
        uint32_t v = png_buf[i] << 16;
        if (i + 1 < png_buf.size()) v |= png_buf[i+1] << 8;
        out += b64[(v >> 18) & 63];
        out += b64[(v >> 12) & 63];
        out += (i + 1 < png_buf.size()) ? b64[(v >> 6) & 63] : '=';
        out += '=';
    }
    return out;
}

// Rasterize a scatter plot into an RGB pixel buffer.
// Each molecule is one pixel; order matters (last drawn wins on collision).
// width_px: desired output width. Height is set to preserve aspect ratio.
static RasterViewport make_viewport(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int width_px
) {
    double xmin = *std::min_element(x.begin(), x.end());
    double xmax = *std::max_element(x.begin(), x.end());
    double ymin = *std::min_element(y.begin(), y.end());
    double ymax = *std::max_element(y.begin(), y.end());
    double xrange = std::max(xmax - xmin, 1.0);
    double yrange = std::max(ymax - ymin, 1.0);

    double margin = 0.01;
    xmin -= xrange * margin; xmax += xrange * margin; xrange *= 1 + 2 * margin;
    ymin -= yrange * margin; ymax += yrange * margin; yrange *= 1 + 2 * margin;

    int height_px = std::max(1, static_cast<int>(width_px * yrange / xrange));
    if (height_px > width_px * 4) height_px = width_px * 4;

    return RasterViewport{xmin, xmax, ymin, ymax, xrange, yrange, width_px, height_px};
}

static inline std::pair<int, int> map_to_pixel(double x, double y, const RasterViewport& vp) {
    int px = static_cast<int>((x - vp.xmin) / vp.xrange * vp.width_px);
    int py = static_cast<int>((vp.ymax - y) / vp.yrange * vp.height_px);
    px = std::max(0, std::min(vp.width_px - 1, px));
    py = std::max(0, std::min(vp.height_px - 1, py));
    return {px, py};
}

static int auto_point_radius_px(const RasterViewport& vp, int n_points) {
    if (n_points <= 0) return 1;
    double area_per_point = static_cast<double>(vp.width_px) * vp.height_px / static_cast<double>(n_points);
    double radius = 0.23 * std::sqrt(area_per_point / PI);
    return std::max(1, std::min(4, static_cast<int>(std::round(radius))));
}

static void draw_disc(
    std::vector<uint8_t>& pixels,
    int width_px,
    int height_px,
    int cx,
    int cy,
    int radius_px,
    uint8_t r, uint8_t g, uint8_t b
) {
    int r2 = radius_px * radius_px;
    for (int dy = -radius_px; dy <= radius_px; ++dy) {
        int yy = cy + dy;
        if (yy < 0 || yy >= height_px) continue;
        for (int dx = -radius_px; dx <= radius_px; ++dx) {
            if (dx * dx + dy * dy > r2) continue;
            int xx = cx + dx;
            if (xx < 0 || xx >= width_px) continue;
            int off = (yy * width_px + xx) * 3;
            pixels[off] = r;
            pixels[off + 1] = g;
            pixels[off + 2] = b;
        }
    }
}

static void draw_line(
    std::vector<uint8_t>& pixels,
    int width_px,
    int height_px,
    int x0, int y0, int x1, int y1,
    uint8_t r, uint8_t g, uint8_t b
) {
    int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    while (true) {
        if (x0 >= 0 && x0 < width_px && y0 >= 0 && y0 < height_px) {
            int off = (y0 * width_px + x0) * 3;
            pixels[off] = r;
            pixels[off + 1] = g;
            pixels[off + 2] = b;
        }
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

static void overlay_polygons(
    std::vector<uint8_t>& pixels,
    const RasterViewport& vp,
    const PolygonCollection& polygons
) {
    for (const auto& [_, poly] : polygons) {
        int n_vertices = 0;
        bool cols_are_vertices = false;
        if (poly.rows() == 2 && poly.cols() >= 2) {
            n_vertices = static_cast<int>(poly.cols());
            cols_are_vertices = true;
        } else if (poly.cols() == 2 && poly.rows() >= 2) {
            n_vertices = static_cast<int>(poly.rows());
            cols_are_vertices = false;
        } else {
            continue;
        }

        for (int i = 0; i < n_vertices; ++i) {
            int j = (i + 1) % n_vertices;
            double x0 = cols_are_vertices ? poly(0, i) : poly(i, 0);
            double y0 = cols_are_vertices ? poly(1, i) : poly(i, 1);
            double x1 = cols_are_vertices ? poly(0, j) : poly(j, 0);
            double y1 = cols_are_vertices ? poly(1, j) : poly(j, 1);
            auto [px0, py0] = map_to_pixel(x0, y0, vp);
            auto [px1, py1] = map_to_pixel(x1, y1, vp);
            draw_line(pixels, vp.width_px, vp.height_px, px0, py0, px1, py1, 0, 0, 0);
        }
    }
}

static std::string render_scatter_impl(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int n,
    int width_px,
    const PolygonCollection* polygons,
    int point_radius_px,
    std::function<void(int, uint8_t&, uint8_t&, uint8_t&)> color_fn
) {
    if (n == 0) return {};

    auto vp = make_viewport(x, y, width_px);
    int point_radius = (point_radius_px > 0) ? point_radius_px : auto_point_radius_px(vp, n);

    // White background.
    std::vector<uint8_t> pixels(vp.width_px * vp.height_px * 3, 255);

    // Writing different pixels from multiple threads is safe.  Two molecules
    // landing on the exact same pixel produce a benign write race (one color
    // wins non-deterministically), which is visually indistinguishable from the
    // serial "last molecule drawn wins" behaviour.
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        auto [px, py] = map_to_pixel(x[i], y[i], vp);
        uint8_t r, g, b;
        color_fn(i, r, g, b);
        draw_disc(pixels, vp.width_px, vp.height_px, px, py, point_radius, r, g, b);
    }

    if (polygons && !polygons->empty()) {
        overlay_polygons(pixels, vp, *polygons);
    }

    return pixels_to_base64_png(pixels, vp.width_px, vp.height_px);
}

std::string render_scatter_png(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::string>& colors,
    const PolygonCollection* polygons,
    int width_px,
    int point_radius_px
) {
    int n = static_cast<int>(x.size());
    // Pre-parse hex colors to avoid repeated string parsing in the hot loop.
    std::vector<uint8_t> cr(n), cg(n), cb(n);
    for (int i = 0; i < n; ++i) hex_to_rgb(colors[i], cr[i], cg[i], cb[i]);

    return render_scatter_impl(x, y, n, width_px, polygons, point_radius_px,
        [&](int i, uint8_t& r, uint8_t& g, uint8_t& b) {
            r = cr[i]; g = cg[i]; b = cb[i];
        });
}

// Blue-orange colormap matching Vega-Lite's "blueorange" scheme.
// confidence=0 → deep orange (#8c3200), confidence=1 → deep blue (#004481)
// We use a simple 3-stop linear interpolation through white at 0.5.
static void blueorange_color(double conf,
                              uint8_t& r, uint8_t& g, uint8_t& b) {
    conf = std::max(0.0, std::min(1.0, conf));
    if (conf < 0.5) {
        // orange end: #8c3200 at 0 → white #ffffff at 0.5
        double t = conf * 2.0;
        r = static_cast<uint8_t>(0x8c + t * (0xff - 0x8c));
        g = static_cast<uint8_t>(0x32 + t * (0xff - 0x32));
        b = static_cast<uint8_t>(0x00 + t * (0xff - 0x00));
    } else {
        // blue end: white #ffffff at 0.5 → #004481 at 1
        double t = (conf - 0.5) * 2.0;
        r = static_cast<uint8_t>(0xff + t * (0x00 - 0xff));
        g = static_cast<uint8_t>(0xff + t * (0x44 - 0xff));
        b = static_cast<uint8_t>(0xff + t * (0x81 - 0xff));
    }
}

std::string render_confidence_png(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& confidence,
    int width_px,
    int point_radius_px
) {
    int n = static_cast<int>(x.size());
    return render_scatter_impl(x, y, n, width_px, nullptr, point_radius_px,
        [&](int i, uint8_t& r, uint8_t& g, uint8_t& b) {
            blueorange_color(confidence[i], r, g, b);
        });
}

// ============================================================================
// vega_noise_histogram
// ============================================================================

nlohmann::json vega_noise_histogram(
    const std::vector<double>& edge_lengths,
    const std::vector<double>& confidence,
    double signal_mu, double signal_sigma,
    double noise_mu, double noise_sigma,
    int nn_id,
    int n_bins
) {
    int n = static_cast<int>(edge_lengths.size());
    if (n == 0) return {};

    // Compute 99th percentile as x_max
    std::vector<double> sorted_el = edge_lengths;
    std::sort(sorted_el.begin(), sorted_el.end());
    double x_max = sorted_el[static_cast<int>(0.99 * (n - 1))];

    // Component weights
    double n1 = 0.0;
    for (double c : confidence) n1 += c;
    double n2 = n - n1;
    double w1 = n1 / (n1 + n2);
    double w2 = n2 / (n1 + n2);

    // Build histogram
    double bin_width = x_max / n_bins;
    std::vector<int> hist(n_bins, 0);
    int n_in_range = 0;
    for (double el : edge_lengths) {
        if (el < x_max) {
            int bin = std::min(static_cast<int>(el / bin_width), n_bins - 1);
            hist[bin]++;
            n_in_range++;
        }
    }

    // Build data: histogram bars + fitted PDF lines
    nlohmann::json values = nlohmann::json::array();
    for (int i = 0; i < n_bins; ++i) {
        double x_center = (i + 0.5) * bin_width;
        double density = static_cast<double>(hist[i]) / (n_in_range * bin_width);

        values.push_back({
            {"x", x_center},
            {"density", density},
            {"type", "Observed"}
        });
        values.push_back({
            {"x", x_center},
            {"density", w1 * normal_pdf(x_center, signal_mu, signal_sigma)}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
            {"type", "Intracellular"}
        });
        values.push_back({
            {"x", x_center},
            {"density", w2 * normal_pdf(x_center, noise_mu, noise_sigma)}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
            {"type", "Background"}
        });
    }

    std::string x_title = "Distance to " + std::to_string(nn_id) + "th nearest neighbor";

    // Color scale shared across all layers — same domain/range in each so that
    // Vega-Lite merges them into a single legend via resolve.
    nlohmann::json color_scale = {
        {"domain", {"Observed", "Intracellular", "Background"}},
        {"range",  {"#1f77b4", "#2ca02c",        "#d62728"  }}
    };
    // Bar layer gets the legend; line layers suppress it (resolve merges all into one).
    nlohmann::json color_with_legend = {
        {"field", "type"}, {"type", "nominal"}, {"scale", color_scale},
        {"legend", {{"title", "Component"}}}
    };
    nlohmann::json color_no_legend = {
        {"field", "type"}, {"type", "nominal"}, {"scale", color_scale},
        {"legend", nullptr}
    };

    return {
        {"$schema", "https://vega.github.io/schema/vega-lite/v5.json"},
        {"title", "Noise estimation"},
        {"width", 500}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"height", 300}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"data", {{"values", values}}},
        // Merge colour legends from all layers into one.
        {"resolve", {{"legend", {{"color", "shared"}}}}},
        {"layer", {
            // Histogram bars — legend shown here
            {
                {"transform", {{{"filter", "datum.type == 'Observed'"}}}},
                {"mark", {{"type", "bar"}, {"opacity", 0.5}}}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
                {"encoding", {
                    {"x", {{"field", "x"}, {"type", "quantitative"}, {"title", x_title},
                           {"bin", {{"binned", true}, {"step", bin_width}}}}},
                    {"x2", {{"expr", "datum.x + " + std::to_string(bin_width)}}},
                    {"y", {{"field", "density"}, {"type", "quantitative"}, {"title", "Density"}}},
                    {"color", color_with_legend}
                }}
            },
            // Signal PDF line
            {
                {"transform", {{{"filter", "datum.type == 'Intracellular'"}}}},
                {"mark", {{"type", "line"}, {"strokeWidth", 3}}}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
                {"encoding", {
                    {"x", {{"field", "x"}, {"type", "quantitative"}}},
                    {"y", {{"field", "density"}, {"type", "quantitative"}}},
                    {"color", color_no_legend}
                }}
            },
            // Noise PDF line
            {
                {"transform", {{{"filter", "datum.type == 'Background'"}}}},
                {"mark", {{"type", "line"}, {"strokeWidth", 3}}}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
                {"encoding", {
                    {"x", {{"field", "x"}, {"type", "quantitative"}}},
                    {"y", {{"field", "density"}, {"type", "quantitative"}}},
                    {"color", color_no_legend}
                }}
            }
        }}
    };
}

// ============================================================================
// vega_gene_frequency
// ============================================================================

nlohmann::json vega_gene_frequency(
    const std::vector<int>& genes,
    const std::vector<double>& confidence,
    const std::vector<std::string>& gene_names
) {
    int n_genes = static_cast<int>(gene_names.size());
    int n = static_cast<int>(genes.size());

    // Count real vs noise per gene
    std::vector<int> real_counts(n_genes, 0);
    std::vector<int> noise_counts(n_genes, 0);
    for (int i = 0; i < n; ++i) {
        int g = genes[i] - 1; // 1-based to 0-based
        if (g < 0 || g >= n_genes) continue;
        if (confidence[i] >= 0.5) {
            real_counts[g]++;
        } else {
            noise_counts[g]++;
        }
    }

    // Filter to genes with >1% of total
    int total = n;
    double threshold = 0.01 * total;

    // Sort genes alphabetically
    std::vector<int> order(n_genes);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return gene_names[a] < gene_names[b];
    });

    nlohmann::json values = nlohmann::json::array();
    for (int g : order) {
        if (real_counts[g] + noise_counts[g] < threshold) continue;
        values.push_back({{"gene", gene_names[g]}, {"count", real_counts[g]}, {"type", "Real"}});
        values.push_back({{"gene", gene_names[g]}, {"count", noise_counts[g]}, {"type", "Noise"}});
    }

    return {
        {"$schema", "https://vega.github.io/schema/vega-lite/v5.json"},
        {"title", "Gene frequency"},
        {"width", 600}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"height", 300}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"data", {{"values", values}}},
        {"mark", "bar"},
        {"encoding", {
            {"x", {{"field", "gene"}, {"type", "nominal"}, {"sort", nullptr},
                   {"axis", {{"labelAngle", -45}}}, {"title", "Gene"}}}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
            {"y", {{"field", "count"}, {"type", "quantitative"}, {"title", "Num. molecules"}}},
            {"color", {
                {"field", "type"}, {"type", "nominal"},
                {"scale", {{"domain", {"Real", "Noise"}}, {"range", {"#1f77b4", "#ff7f0e"}}}},
                {"legend", {{"title", "Type"}}}
            }},
            {"order", {{"field", "type"}}}
        }}
    };
}

// ============================================================================
// vega_gene_structure
// ============================================================================

nlohmann::json vega_gene_structure(const GeneStructureEmbedding& emb) {
    int n = static_cast<int>(emb.x.size());
    nlohmann::json values = nlohmann::json::array();
    for (int i = 0; i < n; ++i) {
        values.push_back({
            {"x", emb.x[i]},
            {"y", emb.y[i]},
            {"gene", emb.gene_names[i]},
            {"size", std::max(emb.marker_sizes[i], 1.0)} // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        });
    }

    // Shared axis encodings: zero:false so the UMAP range isn't forced to
    // include 0 (which would bunch all points in one corner).
    nlohmann::json x_enc = {
        {"field", "x"}, {"type", "quantitative"},
        {"axis", nullptr},
        {"scale", {{"zero", false}}}
    };
    nlohmann::json y_enc = {
        {"field", "y"}, {"type", "quantitative"},
        {"axis", nullptr},
        {"scale", {{"zero", false}}}
    };

    return {
        {"$schema", "https://vega.github.io/schema/vega-lite/v5.json"},
        {"title", "Gene structure"},
        {"width", 500}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"height", 500}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"data", {{"values", values}}},
        {"layer", {
            // Dots
            {
                {"mark", {{"type", "point"}, {"filled", true}, {"opacity", 0.8}}}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
                {"encoding", {
                    {"x", x_enc},
                    {"y", y_enc},
                    {"size", {{"field", "size"}, {"type", "quantitative"},
                              {"scale", {{"range", {20, 400}}}}, {"legend", nullptr}}}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
                    {"tooltip", {
                        {{"field", "gene"}, {"type", "nominal"}},
                        {{"field", "size"}, {"type", "quantitative"}, {"title", "log(count)"}}
                    }}
                }}
            },
            // Gene name labels
            {
                {"mark", {{"type", "text"}, {"dy", -9}, {"fontSize", 10}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
                          {"fontWeight", "normal"}}},
                {"encoding", {
                    {"x", x_enc},
                    {"y", y_enc},
                    {"text", {{"field", "gene"}, {"type", "nominal"}}}
                }}
            }
        }}
    };
}

// ============================================================================
// generate_preview_html
// ============================================================================

std::string generate_preview_html(
    const MoleculeData& data,
    const std::vector<std::string>& gene_colors,
    const std::vector<double>& edge_lengths,
    const NoiseFitResult& noise_result,
    int confidence_nn_id,
    const GeneStructureEmbedding* gene_structure
) {
    // Render PNG images (can be slow — done before HTML assembly)
    std::string scatter_png = render_scatter_png(data.x, data.y, gene_colors);
    std::string conf_png    = render_confidence_png(data.x, data.y, data.confidence);

    // Generate Vega-Lite specs for smaller charts
    auto noise_spec = vega_noise_histogram(
        edge_lengths, data.confidence,
        noise_result.signal_mu, noise_result.signal_sigma,
        noise_result.noise_mu, noise_result.noise_sigma,
        confidence_nn_id
    );
    auto freq_spec = vega_gene_frequency(data.gene, data.confidence, data.gene_names);

    // Compute noise statistics
    int n = data.n_molecules();
    int n_minimal_noise = 0;
    double noise_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        if (data.confidence[i] < 0.01) n_minimal_noise++;
        noise_sum += (1.0 - data.confidence[i]);
    }
    double minimal_noise_pct = (n > 0) ? 100.0 * n_minimal_noise / n : 0.0;
    double expected_noise_pct = (n > 0) ? 100.0 * noise_sum / n : 0.0;

    std::ostringstream html;
    html << R"(<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Baysor Preview Report</title>
<script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 20px; background: #fff; }
h1 { border-bottom: 2px solid #333; padding-bottom: 5px; }
h2 { color: #333; }
.plot-container { margin: 20px 0; }
.png-plot { display: block; max-width: 100%; height: auto; cursor: zoom-in; border: 1px solid #ddd; }
.stats { background: #f5f5f5; padding: 10px 15px; border-radius: 5px; margin: 10px 0; font-size: 14px; }
ul { line-height: 2; }
a { color: #1f77b4; }
</style>
<script>
function openFullRes(src) {
    var parts = src.split(',');
    var mime = parts[0].match(/:(.*?);/)[1];
    var raw = atob(parts[1]);
    var arr = new Uint8Array(raw.length);
    for (var i = 0; i < raw.length; i++) arr[i] = raw.charCodeAt(i);
    var blob = new Blob([arr], {type: mime});
    window.open(URL.createObjectURL(blob), '_blank');
}
</script>
</head>
<body>

<h2>Contents</h2>
<ul>
<li><a href="#transcript_plots">Transcript plots</a></li>
<li><a href="#noise_level">Noise level</a></li>
<li><a href="#gene_structure">Gene structure</a></li>
</ul>

<hr>
<h1 id="transcript_plots">Transcript plots</h1>
<p>Click image to open full resolution in a new tab.</p>
<div class="plot-container">
<img class="png-plot" id="scatter_img" src=")";
    html << scatter_png;
    html << "\" alt=\"Local expression similarity (NCV)\">\n"
            "</div>\n\n"
            "<hr>\n"
            "<h1 id=\"noise_level\">Noise level</h1>\n"
            "<p>Click image to open full resolution in a new tab.</p>\n"
            "<div class=\"plot-container\">\n"
            "<img class=\"png-plot\" id=\"conf_img\" src=\"";
    html << conf_png;
    html << R"(" alt="Transcript confidence">
</div>
<div class="plot-container" id="vg_noise_dist"></div>
<div class="stats">
)";

    html << "Minimal noise level: " << std::fixed << std::setprecision(1) << minimal_noise_pct << "%. "
         << "Expected noise level: " << expected_noise_pct << "%.";

    html << R"(
</div>

<hr>
<h1 id="gene_structure">Gene structure</h1>
<div class="plot-container" id="vg_gene_freq"></div>
)";

    if (gene_structure && !gene_structure->x.empty()) {
        html << R"(<div class="plot-container" id="vg_gene_structure"></div>
)";
    }

    html << R"(
<script>
)";

    // Wire up click-to-full-res for PNG images
    html << R"(
document.querySelectorAll('.png-plot').forEach(function(img) {
    img.addEventListener('click', function() { openFullRes(this.src); });
});
)";

    html << "vegaEmbed('#vg_noise_dist', " << noise_spec.dump() << ", {mode: 'vega-lite'});\n";
    html << "vegaEmbed('#vg_gene_freq', " << freq_spec.dump() << ", {mode: 'vega-lite'});\n";

    if (gene_structure && !gene_structure->x.empty()) {
        auto gs_spec = vega_gene_structure(*gene_structure);
        html << "vegaEmbed('#vg_gene_structure', " << gs_spec.dump() << ", {mode: 'vega-lite'});\n";
    }

    html << R"(
</script>
</body>
</html>
)";

    return html.str();
}

} // namespace baysor
