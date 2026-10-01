// End-to-end tests of the CLI entry point src/cli/main.cpp. main.cpp is only
// compiled into the `baysor` executable, so these tests run that binary
// (BAYSOR_CLI_PATH) as a subprocess on tiny synthetic datasets and check its
// exit code, output and files.

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <tiffio.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "baysor/utils/thread_pool.h"
#include "test_cov_helpers.h"

#if !defined(BAYSOR_CLI_PATH) || defined(_WIN32)

TEST(Cov5Cli, SubprocessTestsNeedPosixAndCliPath) {
    GTEST_SKIP() << "CLI tests need BAYSOR_CLI_PATH and a POSIX shell";
}

#else

namespace {

namespace fs = std::filesystem;
using baysor_test::TempDir;
using baysor_test::csv_column;
using baysor_test::read_text_file;
using baysor_test::cli::CliResult;
using baysor_test::cli::run_cli;

::testing::AssertionResult cli_success(const CliResult& r) {
    if (r.exit_code == 0) return ::testing::AssertionSuccess();
    return ::testing::AssertionFailure()
           << "expected exit code 0, got " << r.exit_code
           << "\n--- stdout ---\n" << r.out << "\n--- stderr ---\n" << r.err;
}

// ---------------------------------------------------------------------------
// Synthetic datasets
// ---------------------------------------------------------------------------

struct MoleculeTable {
    std::vector<double> x, y, z;  // z is empty for 2D tables
    std::vector<std::string> gene;
    std::vector<std::string> cell;  // prior-segmentation labels

    int size() const { return static_cast<int>(x.size()); }
};

// 4 clumps of `per_clump` molecules with 4 cycling genes.
// - with_z: constant z per clump (0/4/8/12), making the table 3D;
// - n_unassigned_per_clump: that many molecules per clump get prior label "0";
// - wide: clumps in one row along x, so the 6000 px wide HTML renders stay small;
// - collinear_last_clump: the last clump lies on one horizontal line, giving a
//   degenerate cell (NaN area-based stats).
MoleculeTable make_clumped_table(int per_clump, bool with_z, int n_unassigned_per_clump,
                                 bool wide = false, bool collinear_last_clump = false) {
    static const char* kGenes[] = {"GeneA", "GeneB", "GeneC", "GeneD"};
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> jit_x(-2.0, 2.0);
    std::uniform_real_distribution<double> jit_y(wide ? -1.0 : -2.0, wide ? 1.0 : 2.0);

    MoleculeTable t;
    for (int c = 0; c < 4; ++c) {
        const double cx = wide ? (10.0 + 200.0 * c) : (10.0 + 20.0 * (c % 2));
        const double cy = wide ? 10.0 : (10.0 + 20.0 * (c / 2));
        const bool collinear = collinear_last_clump && c == 3;
        for (int i = 0; i < per_clump; ++i) {
            t.gene.push_back(kGenes[t.size() % 4]);
            t.x.push_back(cx + jit_x(rng));
            t.y.push_back(collinear ? cy : (cy + jit_y(rng)));
            if (with_z) t.z.push_back(4.0 * c);
            t.cell.push_back(i < n_unassigned_per_clump ? "0" : ("cell" + std::to_string(c + 1)));
        }
    }
    return t;
}

fs::path write_table_csv(const TempDir& tmp, const std::string& name, const MoleculeTable& t) {
    std::ofstream f(tmp.path / name);
    f << "x,y,gene,cell_id" << (t.z.empty() ? "" : ",z") << "\n";
    for (int i = 0; i < t.size(); ++i) {
        f << t.x[i] << "," << t.y[i] << "," << t.gene[i] << "," << t.cell[i];
        if (!t.z.empty()) f << "," << t.z[i];
        f << "\n";
    }
    return tmp.path / name;
}

fs::path write_2d_csv(const TempDir& tmp, int per_clump, bool wide = false,
                      bool collinear_last = false) {
    return write_table_csv(tmp, "mols.csv", make_clumped_table(per_clump, false, 5, wide, collinear_last));
}

// Wide-aspect table for tests whose pipeline renders HTML PNGs.
fs::path write_wide_csv(const TempDir& tmp, int per_clump) {
    return write_2d_csv(tmp, per_clump, /*wide=*/true);
}

// Binary TIFF mask whose foreground covers every molecule pixel (dilated by
// 1), so each clump becomes one labelled component.
fs::path write_mask_tiff(const fs::path& path, const MoleculeTable& t, uint32_t size = 64) {
    std::vector<uint8_t> px(static_cast<size_t>(size) * size, 0);
    for (int i = 0; i < t.size(); ++i) {
        const int col = static_cast<int>(std::round(t.x[i])) - 1;
        const int row = static_cast<int>(std::round(t.y[i])) - 1;
        for (int rr = std::max(row - 1, 0); rr <= std::min<int>(row + 1, size - 1); ++rr)
            for (int cc = std::max(col - 1, 0); cc <= std::min<int>(col + 1, size - 1); ++cc)
                px[static_cast<size_t>(rr) * size + cc] = 1;
    }

    TIFF* tif = TIFFOpen(path.string().c_str(), "w");
    if (!tif) return path;
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, size);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, size);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, size);
    for (uint32_t row = 0; row < size; ++row) {
        TIFFWriteScanline(tif, px.data() + static_cast<size_t>(row) * size, row, 0);
    }
    TIFFClose(tif);
    return path;
}

void write_transcripts_parquet(const fs::path& path, const MoleculeTable& t) {
    arrow::DoubleBuilder xb, yb;
    arrow::StringBuilder gb;
    for (int i = 0; i < t.size(); ++i) {
        EXPECT_TRUE(xb.Append(t.x[i]).ok());
        EXPECT_TRUE(yb.Append(t.y[i]).ok());
        EXPECT_TRUE(gb.Append(t.gene[i]).ok());
    }
    std::shared_ptr<arrow::Array> xa, ya, ga;
    EXPECT_TRUE(xb.Finish(&xa).ok());
    EXPECT_TRUE(yb.Finish(&ya).ok());
    EXPECT_TRUE(gb.Finish(&ga).ok());
    auto schema = arrow::schema({arrow::field("x", arrow::float64()), arrow::field("y", arrow::float64()),
                                 arrow::field("gene", arrow::utf8())});
    auto table = arrow::Table::Make(schema, {xa, ya, ga});

    auto sink = arrow::io::FileOutputStream::Open(path.string());
    ASSERT_TRUE(sink.ok()) << sink.status().ToString();
    EXPECT_TRUE(parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), *sink, 3).ok());
    EXPECT_TRUE((*sink)->Close().ok());
}

// ---------------------------------------------------------------------------
// Output-file assertions
// ---------------------------------------------------------------------------

void expect_files(const fs::path& dir, const std::vector<std::string>& names) {
    for (const auto& name : names) {
        EXPECT_TRUE(fs::is_regular_file(dir / name) && fs::file_size(dir / name) > 0)
            << "expected non-empty output file: " << (dir / name);
    }
}

std::string first_line(const fs::path& p) {
    std::ifstream f(p);
    std::string line;
    std::getline(f, line);
    return line;
}

// Value of "key = value" in a TOML dump.
std::string toml_value(const std::string& toml, const std::string& key) {
    const std::string prefix = key + " = ";
    size_t pos = toml.find(prefix);
    if (pos == std::string::npos) return "";
    pos += prefix.size();
    return toml.substr(pos, toml.find('\n', pos) - pos);
}

bool contains(const std::vector<std::string>& values, const std::string& v) {
    return std::find(values.begin(), values.end(), v) != values.end();
}

}  // namespace

// ============================================================================
// --help / --version
// ============================================================================

TEST(Cov5CliHelp, HelpListsSubcommandsAndOptions) {
    TempDir tmp("help");
    const struct {
        const char* args;
        std::vector<std::string> expected;
    } cases[] = {
        {"--help", {"Baysor", "Run cell segmentation", "preview", "segfree"}},
        {"run --help", {"--output-style", "--count-matrix-format", "--nuclei-genes"}},
    };
    for (const auto& c : cases) {
        const auto r = run_cli(tmp, c.args);
        EXPECT_EQ(r.exit_code, 0) << c.args;
        for (const auto& text : c.expected) EXPECT_NE(r.out.find(text), std::string::npos) << c.args << ": " << text;
    }
}

TEST(Cov5CliHelp, VersionFlags) {
    TempDir tmp("version");
    auto r = run_cli(tmp, "--version");
    EXPECT_EQ(r.exit_code, 0) << r.err;
    EXPECT_EQ(r.out, std::string("baysor ") + BAYSOR_VERSION + "\n");
    for (const auto* command : {"run", "preview", "segfree"}) {
        r = run_cli(tmp, std::string(command) + " --version");
        EXPECT_EQ(r.exit_code, 0) << command << ": " << r.err;
        // Sopa parses this with packaging.version.Version: no "baysor " prefix.
        EXPECT_EQ(r.out, std::string(BAYSOR_VERSION) + "\n") << command;
    }
}

TEST(Cov5Cli, ThreadCountFromOmpNumThreads) {
    // The thread count is logged before the (missing) input is read.
    TempDir tmp("threads");
    const struct { const char* env; const char* args; const char* expected; } cases[] = {
        {"3", "", "Using 3 threads"},
        {"3", "-t 2 ", "Using 2 threads"},  // --threads wins
        {"5,2", "", "Using 5 threads"},     // first entry of an OpenMP list
    };
    for (const auto& c : cases) {
        setenv("OMP_NUM_THREADS", c.env, 1);
        const auto r = run_cli(tmp, std::string("run ") + c.args + "-s 5 '" +
                                        (tmp.path / "missing.csv").string() + "'");
        EXPECT_NE(r.out.find(c.expected), std::string::npos) << c.env << " " << c.args << r.out;
    }
    setenv("OMP_NUM_THREADS", "bogus", 1);  // ignored: falls back to the core count
    const auto r = run_cli(tmp, "run -s 5 '" + (tmp.path / "missing.csv").string() + "'");
    EXPECT_NE(r.out.find("Using " + std::to_string(baysor::default_thread_count()) + " threads"),
              std::string::npos) << r.out;
    unsetenv("OMP_NUM_THREADS");
}

// ============================================================================
// Invalid invocations: CLI11 parse errors, option validation, input errors
// and guarded paths of cmd_run / cmd_preview
// ============================================================================

TEST(Cov5Cli, InvalidInvocationsFailWithMessage) {
    TempDir tmp("cli_errors");
    const std::string csv = "'" + write_2d_csv(tmp, 20).string() + "'";
    const std::string wide_csv = "'" + write_table_csv(tmp, "wide.csv", make_clumped_table(20, false, 5, true)).string() + "'";
    const std::string dir = tmp.path.string();
    // Boundary polygon that overlaps the molecule bounding box but contains no
    // molecules (the clumps occupy [8,12] and [28,32] on both axes).
    tmp.write("bounds.csv", "cell_id,vertex_x,vertex_y\ncell1,16,16\ncell1,16,24\ncell1,24,24\ncell1,24,16\n");
    tmp.write("blocked", "i am a file\n");
    fs::create_directories(tmp.path / "xen");
    tmp.write("xen/experiment.xenium", "run_folder,panel\nsynthetic,none\n");

    const struct {
        std::string args;
        int exit_code;  // -1: any non-zero code (CLI11 parse errors)
        std::vector<std::string> messages;  // expected in stdout + stderr
    } cases[] = {
        {"", -1, {"A subcommand is required"}},
        {"run " + csv + " -m 10 -s 2.5 --frobnicate", -1, {"--frobnicate"}},
        {"run -m 10 -s 2.5", -1, {"coordinates"}},
        {"run " + csv + " -m 10 -s 2.5 --iters abc", -1, {"iters"}},
        {"run " + csv + " -m 10 -s 2.5 --polygon-format geo -o '" + dir + "/seg_poly'", -1, {"polygon-format"}},
        {"run " + csv + " -m 10 -s 2.5 --output-style xml", 1, {"Unknown output style: xml"}},
        {"run " + csv + " -m 10 -s 2.5 --cluster-method kmeans", 1, {"cluster_method must be one of"}},
        {"run " + csv + " -m 10", 1, {"Either prior_segmentation or --scale must be provided."}},
        {"run " + csv + " -c '" + dir + "/missing.toml' -m 10 -s 2.5", 1, {"Failed to load config"}},
        {"run '" + dir + "/nope.csv' -m 10 -s 2.5", 1,
         {"Arrow error: IOError: Failed to open local file", "No such file or directory"}},
        {"run " + csv + " '" + dir + "/missing_mask.tiff' -m 10 -s 2.5", 1, {"Cannot open TIFF mask file"}},
        {"run '" + dir + "/xen/experiment.xenium' -m 10 -s 2.5", 1, {"Could not locate transcripts"}},
        // The boundary contains no molecules: scale estimation warns and the
        // run aborts before segmentation.
        {"run " + csv + " '" + dir + "/bounds.csv' -m 10 -o '" + dir + "/seg'", 1,
         {"Could not estimate scale from prior", "Scale could not be determined"}},
        {"run " + csv + " -m 10 -s 2.5 --nuclei-genes GeneA -o '" + dir + "/seg'", 1, {"not yet implemented"}},
        {"run " + csv + " -m 10 -s 2.5 -o '" + dir + "/blocked/seg'", 1, {"Could not create output directory"}},
        {"preview " + wide_csv + " -m 10 -o '" + dir + "/no_such_dir/preview.html'", 1, {"Could not write to"}},
    };
    for (const auto& c : cases) {
        const auto r = run_cli(tmp, c.args);
        if (c.exit_code < 0) EXPECT_NE(r.exit_code, 0) << c.args;
        else EXPECT_EQ(r.exit_code, c.exit_code) << c.args;
        for (const auto& msg : c.messages) {
            EXPECT_NE((r.out + r.err).find(msg), std::string::npos) << c.args << "\n" << r.out << r.err;
        }
    }
    EXPECT_FALSE(fs::exists(tmp.path / "seg_poly"));
}

// ============================================================================
// Successful `run` invocations covering the output-style / prior / cluster /
// dimensionality matrix
// ============================================================================

TEST(Cov5CliRun, PriorColumnLegacy2DBundle) {
    TempDir tmp("run_prior_col");
    // The last clump yields a zero-area cell: NaN density/elongation.
    auto csv = write_2d_csv(tmp, 50, /*wide=*/false, /*collinear_last=*/true);
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run '" + csv.string() + "' ':cell_id' -m 10 --iters 12 -o '" + out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    EXPECT_NE(r.out.find("prior-aware n_cells_init="), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("Segmentation complete"), std::string::npos);

    expect_files(out, {"segmentation.csv", "segmentation_cell_stats.csv", "segmentation_counts.loom",
                       "segmentation_polygons_2d.json", "segmentation_params.dump.toml", "segmentation_log.log"});

    // One row per molecule, with the columns that need clustering, NCV
    // colours and assignment confidence.
    EXPECT_EQ(csv_column(out / "segmentation.csv", "cell").size(), 200u);
    const std::string seg_header = first_line(out / "segmentation.csv");
    for (const auto* col : {"ncv_color", "assignment_confidence", "cluster"}) {
        EXPECT_NE(seg_header.find(col), std::string::npos) << seg_header;
    }
    EXPECT_EQ(seg_header.find(",z"), std::string::npos) << seg_header;

    const fs::path stats = out / "segmentation_cell_stats.csv";
    const std::string stats_header = first_line(stats);
    for (const auto* col : {"cluster", "max_cluster_frac", "lifespan", "avg_assignment_confidence"}) {
        EXPECT_NE(stats_header.find(col), std::string::npos) << stats_header;
    }
    EXPECT_GE(csv_column(stats, "density").size(), 2u);
    EXPECT_TRUE(contains(csv_column(stats, "density"), "nan")) << read_text_file(stats);
    EXPECT_TRUE(contains(csv_column(stats, "elongation"), "nan")) << read_text_file(stats);

    EXPECT_NE(read_text_file(out / "segmentation_polygons_2d.json").find("FeatureCollection"), std::string::npos);

    // Params dump records the prior and a positive scale estimated from it.
    const std::string dump = read_text_file(out / "segmentation_params.dump.toml");
    EXPECT_EQ(toml_value(dump, "type"), "\"column\"") << dump;
    EXPECT_GT(std::stod(toml_value(dump, "scale")), 0.0) << dump;
    EXPECT_EQ(toml_value(dump, "iters"), "12") << dump;

    // Dual logger: messages also land in the log file.
    EXPECT_NE(read_text_file(out / "segmentation_log.log").find("Segmentation complete"), std::string::npos);
}

TEST(Cov5CliRun, PriorColumnLegacyGeometryCollectionLegacyHasIntegerCellIds) {
    TempDir tmp("run_legacy_poly");
    // The collinear clump's polygon estimation fails and the fallback keeps
    // it in the polygons file (kharchenkolab/Baysor#165).
    auto csv = write_2d_csv(tmp, 50, /*wide=*/false, /*collinear_last=*/true);
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run '" + csv.string() + "' ':cell_id' -m 10 --iters 12 -s 2.5"
                          " --polygon-format GeometryCollectionLegacy -o '" + out.string() + "'");
    ASSERT_TRUE(cli_success(r));

    const auto doc = nlohmann::json::parse(read_text_file(out / "segmentation_polygons_2d.json"));
    EXPECT_EQ(doc.at("type"), "GeometryCollection");
    std::set<int> poly_ids;
    for (const auto& geom : doc.at("geometries")) {
        ASSERT_TRUE(geom.at("cell").is_number_integer()) << geom.dump();
        poly_ids.insert(geom.at("cell").get<int>());
    }
    // The integer ids are exactly the CSV's cell_<n> names.
    std::set<int> csv_ids;
    for (const auto& cell : csv_column(out / "segmentation.csv", "cell")) {
        if (cell.rfind("cell_", 0) == 0) csv_ids.insert(std::stoi(cell.substr(5)));
    }
    EXPECT_EQ(poly_ids, csv_ids);
    EXPECT_FALSE(poly_ids.empty());
}

TEST(Cov5CliRun, ParquetStyleWithPlotAndIgnoredFormatWarnings) {
    TempDir tmp("run_parquet_plot");
    auto csv = write_wide_csv(tmp, 40);
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run '" + csv.string() + "' -m 10 -s 2.5 --iters 12 -p --output-style parquet"
                          " --cluster-method none --polygon-format GeometryCollection"
                          " --count-matrix-format tsv -o '" + out.string() + "'");
    ASSERT_TRUE(cli_success(r));

    // Non-default polygon/count formats are ignored under parquet output.
    EXPECT_NE(r.out.find("--polygon-format is ignored for output style 'parquet'"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("--count-matrix-format is ignored for output style 'parquet'"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("Generating HTML run report"), std::string::npos);

    expect_files(out, {"molecules.parquet", "cells.parquet", "cell_boundaries.parquet", "feature_matrix.h5",
                       "diagnostic_report.html", "segmentation_plot.html", "run_params.toml", "run.log"});
    EXPECT_NE(read_text_file(out / "diagnostic_report.html").find("html"), std::string::npos);
    EXPECT_NE(read_text_file(out / "segmentation_plot.html").find("html"), std::string::npos);
    EXPECT_FALSE(fs::exists(out / "segmentation.csv"));
    EXPECT_FALSE(fs::exists(out / "segmentation_counts.loom"));
}

TEST(Cov5CliRun, ConfigFileWithCliOverride3DAndLouvain) {
    TempDir tmp("run_config3d");
    auto csv = write_table_csv(tmp, "mols3d.csv", make_clumped_table(30, /*with_z=*/true, 5));
    const auto cfg = tmp.write("baysor.toml",
                               "[molecules]\nmin_molecules_per_cell = 10\n\n"
                               "[segmentation]\nscale = 6.0\ncluster_method = \"louvain\"\niters = 400\n");
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run -c '" + cfg + "' '" + csv.string() +
                          "' --iters 12 --polygon-format GeometryCollection -o '" + out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    EXPECT_NE(r.out.find("Loaded config from"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("Louvain on NCV kNN graph"), std::string::npos) << r.out;

    expect_files(out, {"segmentation.csv", "segmentation_cell_stats.csv", "segmentation_counts.loom",
                       "segmentation_polygons_3d.json", "segmentation_params.dump.toml"});
    // 3D input produces a z column in the cell stats.
    EXPECT_NE(first_line(out / "segmentation_cell_stats.csv").find(",z"), std::string::npos);

    // Config supplies defaults; the explicit --iters flag wins.
    const std::string dump = read_text_file(out / "segmentation_params.dump.toml");
    EXPECT_EQ(toml_value(dump, "cluster_method"), "\"louvain\"") << dump;
    EXPECT_EQ(toml_value(dump, "min_molecules_per_cell"), "10") << dump;
    EXPECT_EQ(toml_value(dump, "iters"), "12") << dump;
    EXPECT_EQ(std::stod(toml_value(dump, "scale")), 6.0) << dump;
}

TEST(Cov5CliRun, Leiden3DParquetBundle) {
    TempDir tmp("run_leiden3d");
    auto csv = write_table_csv(tmp, "mols3d.csv", make_clumped_table(30, /*with_z=*/true, 5));
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run '" + csv.string() + "' -m 10 -s 2.5 --iters 12"
                          " --cluster-method leiden --output-style parquet -o '" + out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    EXPECT_NE(r.out.find("Leiden on NCV kNN graph"), std::string::npos) << r.out;
    // 3D input writes both the combined 2D projection and the per-layer stack.
    expect_files(out, {"molecules.parquet", "cells.parquet", "feature_matrix.h5", "run_params.toml", "run.log",
                       "cell_boundaries.parquet", "cell_boundaries_3d.parquet"});
    const std::string dump = read_text_file(out / "run_params.toml");
    EXPECT_EQ(toml_value(dump, "cluster_method"), "\"leiden\"") << dump;
}

TEST(Cov5CliRun, PlotHtmlWriteFailuresAreReported) {
    TempDir tmp("run_plot_fail");
    auto csv = write_wide_csv(tmp, 20);
    const fs::path out = tmp.path / "seg";
    fs::create_directories(out / "diagnostic_report.html");
    fs::create_directories(out / "segmentation_plot.html");

    auto r = run_cli(tmp, "run '" + csv.string() + "' -m 10 -s 2.5 --iters 10 -p"
                          " --cluster-method none --skip-ncv-color -o '" + out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    EXPECT_NE(r.out.find("Could not write diagnostic report"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("Could not write segmentation plot"), std::string::npos) << r.out;
}

TEST(Cov5CliRun, ImagePriorNoClusteringTsvWithoutPolygons) {
    TempDir tmp("run_image");
    auto table = make_clumped_table(50, /*with_z=*/false, 0);
    auto csv = write_table_csv(tmp, "mols.csv", table);
    auto mask = write_mask_tiff(tmp.path / "mask.tif", table);
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run '" + csv.string() + "' '" + mask.string() +
                          "' -m 10 --iters 12 --cluster-method none --count-matrix-format tsv"
                          " --polygon-format none --skip-ncv-color -o '" + out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    expect_files(out, {"segmentation.csv", "segmentation_cell_stats.csv", "segmentation_counts.tsv",
                       "segmentation_params.dump.toml"});

    // tsv count matrix: header row keyed by gene, one row per gene.
    const std::string tsv = read_text_file(out / "segmentation_counts.tsv");
    EXPECT_EQ(tsv.find("gene\tcell_1"), 0u) << tsv.substr(0, 80);
    EXPECT_NE(tsv.find("GeneA"), std::string::npos) << tsv;

    EXPECT_FALSE(fs::exists(out / "segmentation_polygons_2d.json"));

    // --cluster-method none: no cluster columns; --skip-ncv-color: no ncv_color.
    const std::string stats_header = first_line(out / "segmentation_cell_stats.csv");
    EXPECT_NE(stats_header.find("n_transcripts"), std::string::npos) << stats_header;
    EXPECT_EQ(stats_header.find("max_cluster_frac"), std::string::npos) << stats_header;
    EXPECT_EQ(stats_header.find(",cluster"), std::string::npos) << stats_header;
    EXPECT_EQ(first_line(out / "segmentation.csv").find("ncv_color"), std::string::npos);

    const std::string dump = read_text_file(out / "segmentation_params.dump.toml");
    EXPECT_EQ(toml_value(dump, "type"), "\"image\"") << dump;
    EXPECT_GT(std::stod(toml_value(dump, "scale")), 0.0) << dump;
}

// ============================================================================
// preview / segfree
// ============================================================================

TEST(Cov5CliPreview, WritesHtmlReport) {
    TempDir tmp("preview_ok");
    auto csv = write_wide_csv(tmp, 40);
    const fs::path out = tmp.path / "preview.html";

    auto r = run_cli(tmp, "preview '" + csv.string() + "' -m 10 -o '" + out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    EXPECT_NE(r.out.find("Preview saved to"), std::string::npos) << r.out;
    EXPECT_NE(read_text_file(out).find("html"), std::string::npos);
}

TEST(Cov5CliSegfree, InfersDefaultKAndWritesLoom) {
    TempDir tmp("segfree_default");
    auto csv = write_2d_csv(tmp, 50);
    const fs::path out = tmp.path / "ncvs.loom";

    auto r = run_cli(tmp, "segfree '" + csv.string() + "' -m 10 -o '" + out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    // default_param_value("composition_neighborhood", min_molecules_per_cell=10, n_genes=4) == 10
    EXPECT_NE(r.out.find("Using k=10 neighbors for NCV composition"), std::string::npos) << r.out;
    expect_files(tmp.path, {"ncvs.loom"});
}

TEST(Cov5CliSegfree, ExplicitKWithXeniumManifestInput) {
    TempDir tmp("segfree_xenium");
    fs::create_directories(tmp.path / "xen");
    tmp.write("xen/experiment.xenium", "run_folder,panel\nsynthetic,none\n");
    write_transcripts_parquet(tmp.path / "xen" / "transcripts.parquet", make_clumped_table(50, false, 0));
    const fs::path out = tmp.path / "ncvs.loom";

    auto r = run_cli(tmp, "segfree '" + tmp.file("xen/experiment.xenium") + "' -m 10 -k 12 -o '" + out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    EXPECT_NE(r.out.find("Using k=12 neighbors for NCV composition"), std::string::npos) << r.out;
    expect_files(tmp.path, {"ncvs.loom"});
}

#endif  // BAYSOR_CLI_PATH && !defined(_WIN32)
