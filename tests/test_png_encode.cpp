// PNG images of the HTML reports: decode the base64 data URIs produced by the
// renderers with an independent minimal PNG reader (zlib inflate + all five
// filter types) and check the pixels.

#include "baysor/reporting/preview_report.h"
#include "baysor/reporting/run_report.h"
#include "baysor/utils/options.h"
#include "baysor/utils/thread_pool.h"

#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <zlib.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct DecodedPng {
    uint32_t width = 0, height = 0;
    std::vector<uint8_t> rgb;  // row-major, 3 bytes per pixel
    int n_idat = 0;
};

std::vector<uint8_t> base64_decode(const std::string& s) {
    auto val = [](char c) -> int {
        if (c >= 'A' && c <= 'Z') return c - 'A';
        if (c >= 'a' && c <= 'z') return c - 'a' + 26;
        if (c >= '0' && c <= '9') return c - '0' + 52;
        if (c == '+') return 62;
        if (c == '/') return 63;
        return -1;
    };
    std::vector<uint8_t> out;
    uint32_t acc = 0;
    int bits = 0;
    for (char c : s) {
        const int v = val(c);
        if (v < 0) break;  // '=' padding
        acc = (acc << 6) | static_cast<uint32_t>(v);
        bits += 6;
        if (bits >= 8) {
            bits -= 8;
            out.push_back(static_cast<uint8_t>(acc >> bits));
        }
    }
    return out;
}

uint32_t be32(const uint8_t* p) {
    return (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | uint32_t(p[3]);
}

uint8_t paeth(int a, int b, int c) {
    const int p = a + b - c, pa = std::abs(p - a), pb = std::abs(p - b), pc = std::abs(p - c);
    return static_cast<uint8_t>((pa <= pb && pa <= pc) ? a : (pb <= pc ? b : c));
}

// Minimal reader for 8-bit RGB, non-interlaced PNGs.
DecodedPng decode_png_data_uri(const std::string& uri) {
    const std::string prefix = "data:image/png;base64,";
    if (uri.compare(0, prefix.size(), prefix) != 0) throw std::runtime_error("not a PNG data URI");
    const std::vector<uint8_t> png = base64_decode(uri.substr(prefix.size()));
    static const uint8_t sig[8] = {0x89, 'P', 'N', 'G', '\r', '\n', 0x1a, '\n'};
    if (png.size() < 8 || !std::equal(sig, sig + 8, png.begin())) throw std::runtime_error("bad signature");

    DecodedPng out;
    std::vector<uint8_t> zdata;
    bool seen_iend = false;
    size_t pos = 8;
    while (pos + 12 <= png.size()) {
        const uint32_t len = be32(&png[pos]);
        const std::string type(reinterpret_cast<const char*>(&png[pos + 4]), 4);
        if (pos + 12 + len > png.size()) throw std::runtime_error("truncated chunk");
        const uint8_t* data = &png[pos + 8];
        const uint32_t crc = static_cast<uint32_t>(crc32(crc32(0L, Z_NULL, 0), &png[pos + 4], 4 + len));
        if (crc != be32(&png[pos + 8 + len])) throw std::runtime_error("bad CRC in " + type);
        if (type == "IHDR") {
            out.width = be32(data);
            out.height = be32(data + 4);
            if (data[8] != 8 || data[9] != 2 || data[12] != 0) throw std::runtime_error("unsupported format");
        } else if (type == "IDAT") {
            zdata.insert(zdata.end(), data, data + len);
            ++out.n_idat;
        } else if (type == "IEND") {
            seen_iend = true;
        }
        pos += 12 + len;
    }
    if (!seen_iend || pos != png.size()) throw std::runtime_error("missing IEND or trailing bytes");

    const size_t stride = size_t(out.width) * 3;
    std::vector<uint8_t> raw((stride + 1) * out.height);
    uLongf raw_len = static_cast<uLongf>(raw.size());
    if (uncompress(raw.data(), &raw_len, zdata.data(), static_cast<uLong>(zdata.size())) != Z_OK ||
        raw_len != raw.size()) {
        throw std::runtime_error("inflate failed");
    }

    out.rgb.assign(stride * out.height, 0);
    for (size_t y = 0; y < out.height; ++y) {
        const uint8_t f = raw[y * (stride + 1)];
        const uint8_t* src = &raw[y * (stride + 1) + 1];
        uint8_t* dst = &out.rgb[y * stride];
        const uint8_t* up = y > 0 ? &out.rgb[(y - 1) * stride] : nullptr;
        for (size_t i = 0; i < stride; ++i) {
            const int a = i >= 3 ? dst[i - 3] : 0, b = up ? up[i] : 0, c = (up && i >= 3) ? up[i - 3] : 0;
            int pred = 0;
            switch (f) {
                case 0: pred = 0; break;
                case 1: pred = a; break;
                case 2: pred = b; break;
                case 3: pred = (a + b) / 2; break;
                case 4: pred = paeth(a, b, c); break;
                default: throw std::runtime_error("bad filter type");
            }
            dst[i] = static_cast<uint8_t>(src[i] + pred);
        }
    }
    return out;
}

class PoolSizeGuard {
public:
    explicit PoolSizeGuard(int n) : old_(baysor::thread_pool_size()) { baysor::set_thread_pool_size(n); }
    ~PoolSizeGuard() { baysor::set_thread_pool_size(old_); }
private:
    int old_;
};

baysor::ScatterRaster test_raster(int n_points, int width, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 100.0);
    std::uniform_int_distribution<int> c(0, 15);
    const char* hex = "0123456789abcdef";
    std::vector<double> x(n_points), y(n_points);
    std::vector<std::string> colors(n_points);
    for (int i = 0; i < n_points; ++i) {
        x[i] = u(rng);
        y[i] = u(rng) * 0.7;
        colors[i] = std::string("#") + hex[c(rng)] + hex[c(rng)] + hex[c(rng)] + hex[c(rng)] + hex[c(rng)] + hex[c(rng)];
    }
    baysor::PolygonCollection polygons;
    Eigen::MatrixXd quad(2, 4);
    quad << 10.0, 60.0, 60.0, 10.0,
            10.0, 10.0, 50.0, 50.0;
    polygons["c1"] = quad;
    return baysor::rasterize_scatter(x, y, colors, &polygons, width);
}

} // namespace

TEST(PngEncode, RasterRoundTripsThroughPng) {
    const baysor::ScatterRaster r = test_raster(3000, 517, 1);
    ASSERT_FALSE(r.empty());
    const DecodedPng img = decode_png_data_uri(baysor::encode_png_data_uris({r})[0]);
    EXPECT_EQ(img.width, static_cast<uint32_t>(r.width_px));
    EXPECT_EQ(img.height, static_cast<uint32_t>(r.height_px));
    EXPECT_EQ(img.rgb, r.pixels);
}

TEST(PngEncode, ConcurrentEncodingMatchesSerialAndKeepsOrder) {
    std::vector<baysor::ScatterRaster> rasters = {
        test_raster(2000, 400, 2), baysor::ScatterRaster{}, test_raster(500, 123, 3), test_raster(4000, 700, 4)};
    std::vector<std::string> serial;
    {
        PoolSizeGuard guard(1);
        serial = baysor::encode_png_data_uris(rasters);
    }
    PoolSizeGuard guard(4);
    const std::vector<std::string> parallel = baysor::encode_png_data_uris(rasters);
    ASSERT_EQ(parallel.size(), rasters.size());
    EXPECT_EQ(parallel, serial);
    EXPECT_EQ(parallel[1], "");
    for (size_t i : {0u, 2u, 3u}) {
        EXPECT_EQ(decode_png_data_uri(parallel[i]).rgb, rasters[i].pixels) << "image " << i;
    }
}

TEST(PngEncode, ScatterPngDecodesToTheDrawnColours) {
    // Points on a grid with distinct colours; the raster is mostly white.
    std::vector<double> x, y;
    std::vector<std::string> colors;
    std::set<std::tuple<int, int, int>> expected = {{255, 255, 255}};
    const char* hex = "0123456789abcdef";
    for (int i = 0; i < 40; ++i) {
        x.push_back(i % 8);
        y.push_back(i / 8);
        const int r = (i * 37) % 256, g = (i * 91) % 256, b = (i * 53) % 200;
        colors.push_back(std::string("#") + hex[r >> 4] + hex[r & 15] + hex[g >> 4] + hex[g & 15] +
                         hex[b >> 4] + hex[b & 15]);
        expected.insert({r, g, b});
    }

    const std::string uri = baysor::render_scatter_png(x, y, colors, nullptr, /*width_px=*/333, /*point_radius_px=*/3);
    const DecodedPng img = decode_png_data_uri(uri);
    EXPECT_EQ(img.width, 333u);
    EXPECT_GT(img.height, 100u);
    ASSERT_EQ(img.rgb.size(), size_t(img.width) * img.height * 3);

    std::set<std::tuple<int, int, int>> seen;
    for (size_t p = 0; p < img.rgb.size(); p += 3) seen.insert({img.rgb[p], img.rgb[p + 1], img.rgb[p + 2]});
    EXPECT_EQ(seen, expected);
}

TEST(PngEncode, LargeImageIsSplitIntoSeveralIdatChunks) {
    // Noise-like colours compress poorly, so the deflate stream exceeds the
    // 1 MiB IDAT chunk size and the chunking path is exercised.
    std::vector<double> x, y;
    std::vector<std::string> colors;
    const char* hex = "0123456789abcdef";
    uint32_t state = 12345;
    for (int i = 0; i < 1500 * 1500 / 4; ++i) {
        x.push_back(i % 750);
        y.push_back(i / 750);
        std::string c = "#";
        for (int k = 0; k < 6; ++k) {
            state = state * 1664525u + 1013904223u;
            c += hex[state >> 28];
        }
        colors.push_back(c);
    }
    const DecodedPng img = decode_png_data_uri(
        baysor::render_scatter_png(x, y, colors, nullptr, /*width_px=*/1500, /*point_radius_px=*/1));
    EXPECT_EQ(img.width, 1500u);
    EXPECT_GT(img.n_idat, 1);
    EXPECT_EQ(img.rgb.size(), size_t(img.width) * img.height * 3);
}

TEST(PlotSize, DefaultMatchesPlottingOptions) {
    EXPECT_EQ(baysor::kDefaultMaxPlotSize, baysor::PlottingOptions{}.max_plot_size);
}

TEST(PlotSize, LongerSideIsMaxPlotSize) {
    using baysor::scatter_width_for_max_size;
    auto raster_dims = [](const std::vector<double>& x, const std::vector<double>& y, int max_size) {
        const int w = scatter_width_for_max_size(x, y, max_size);
        const auto r = baysor::rasterize_scatter(x, y, std::vector<std::string>(x.size(), "#000000"), nullptr, w, 1);
        return std::make_pair(r.width_px, r.height_px);
    };
    // Square and wide data: the width is the longer side.
    EXPECT_EQ(raster_dims({0, 10, 0, 10}, {0, 0, 10, 10}, 500), std::make_pair(500, 500));
    auto wide = raster_dims({0, 30}, {0, 10}, 600);
    EXPECT_EQ(wide.first, 600);
    EXPECT_LE(wide.second, 600);
    // Tall data: the height is the longer side and stays within the limit.
    auto tall = raster_dims({0, 10}, {0, 25}, 1000);
    EXPECT_LE(tall.second, 1000);
    EXPECT_GE(tall.second, 995);
    EXPECT_NEAR(tall.first, 400, 2);
    // Taller than 4:1: the height is capped at 4x the width, so the width is max/4.
    auto very_tall = raster_dims({0, 1}, {0, 100}, 1000);
    EXPECT_EQ(very_tall, std::make_pair(250, 1000));
    // Invalid limits fall back to the default; empty data keeps the limit.
    EXPECT_EQ(scatter_width_for_max_size({0, 1}, {0, 1}, 0), baysor::kDefaultMaxPlotSize);
    EXPECT_EQ(scatter_width_for_max_size({}, {}, 700), 700);
}

TEST(PlotSize, SegmentationReportHonoursMaxPlotSize) {
    baysor::MoleculeData d;
    for (int i = 0; i < 200; ++i) {
        d.x.push_back(i % 20);
        d.y.push_back(i / 20 * 4.0);  // 19 x 36: taller than wide
        d.gene.push_back(1 + i % 2);
        d.confidence.push_back(0.9);
    }
    d.gene_names = {"A", "B"};
    std::vector<int> assignment(200, 1);
    std::vector<std::string> ncv(200, "#336699");
    const std::string html = baysor::generate_run_segmentation_html(d, assignment, ncv, nullptr, nullptr, 320);
    const std::regex uri_re("data:image/png;base64,[A-Za-z0-9+/=]+");
    int n_images = 0;
    for (auto it = std::sregex_iterator(html.begin(), html.end(), uri_re); it != std::sregex_iterator(); ++it) {
        const DecodedPng img = decode_png_data_uri(it->str());
        EXPECT_LE(std::max(img.width, img.height), 320u);
        EXPECT_GE(img.height, 315u);
        EXPECT_LT(img.width, img.height);
        ++n_images;
    }
    EXPECT_EQ(n_images, 2);  // assignment and NCV colours
}
