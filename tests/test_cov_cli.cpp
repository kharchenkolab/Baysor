// COV-5: end-to-end tests for the CLI entry point src/cli/main.cpp.
//
// main.cpp is compiled only into the `baysor` executable, so it cannot be
// linked into the unit-test binary. These tests spawn the (instrumented)
// `baysor` binary — whose path CMake injects as BAYSOR_CLI_PATH — as a
// subprocess and assert on its exit code, stdout/stderr and output files.
// When the coverage build runs the suite, the subprocess writes its own
// .gcda counters for src/cli/main.cpp, which the `coverage` target collects.
//
// All datasets are tiny synthetic clumps (3 genes, 4 clumps) and all runs use
// small iteration counts so every CLI invocation finishes quickly.

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <tiffio.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "test_cov_helpers.h"

#ifndef BAYSOR_CLI_PATH

// The CLI path is only injected by the BAYSOR_WITH_TESTS CMake block. Keep a
// placeholder test so the suite stays visible if the define goes missing.
TEST(Cov5Cli, BaysorCliPathAvailable) {
    GTEST_SKIP() << "BAYSOR_CLI_PATH is not defined; CLI end-to-end tests are disabled";
}

#elif defined(_WIN32)

// The subprocess runner below shells out with sh-style quoting and decodes
// exit codes via WEXITSTATUS, so the end-to-end CLI tests are POSIX-only.
TEST(Cov5Cli, SubprocessTestsArePosixOnly) {
    GTEST_SKIP() << "CLI subprocess tests require a POSIX shell and sys/wait.h";
}

#else  // BAYSOR_CLI_PATH && !defined(_WIN32)

#include <sys/wait.h>
#include <unistd.h>

namespace {

namespace fs = std::filesystem;

std::string cli_binary() {
    return BAYSOR_CLI_PATH;
}

// ---------------------------------------------------------------------------
// Per-test unique temp directory (removed on destruction)
// ---------------------------------------------------------------------------

// Portable RAII temp directory (see tests/test_cov_helpers.h): unique via a
// counter plus a random suffix, no getpid()/POSIX.
using TempDir = baysor_test::TempDir;

// ---------------------------------------------------------------------------
// Subprocess runner
// ---------------------------------------------------------------------------

std::string read_text_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

struct CliResult {
    int exit_code = -1;  // -1 = process did not exit normally
    std::string out;     // stdout
    std::string err;     // stderr
};

CliResult run_cli(const TempDir& tmp, const std::string& args) {
    const fs::path out_p = tmp.path / "stdout.txt";
    const fs::path err_p = tmp.path / "stderr.txt";
    const std::string cmd = "'" + cli_binary() + "' " + args +
                            " > '" + out_p.string() + "' 2> '" + err_p.string() + "'";
    const int status = std::system(cmd.c_str());

    CliResult r;
    if (status >= 0 && WIFEXITED(status)) {
        r.exit_code = WEXITSTATUS(status);
    }
    r.out = read_text_file(out_p);
    r.err = read_text_file(err_p);
    return r;
}

::testing::AssertionResult cli_success(const CliResult& r) {
    if (r.exit_code != 0) {
        return ::testing::AssertionFailure()
               << "expected exit code 0, got " << r.exit_code
               << "\n--- stdout ---\n" << r.out
               << "\n--- stderr ---\n" << r.err;
    }
    return ::testing::AssertionSuccess();
}

// ---------------------------------------------------------------------------
// Synthetic datasets
// ---------------------------------------------------------------------------

struct MoleculeTable {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;       // empty for 2D tables
    std::vector<std::string> gene;
    std::vector<std::string> cell;  // prior-segmentation labels

    int size() const { return static_cast<int>(x.size()); }
};

// 4 clumps, 4 genes cycling, jitter +-2 in x ( +-1 in wide mode ).
// NOTE: this table predates the BUG-2 fix; the original 4-gene choice
// (rather than 3) was made because fewer genes than the default mrf
// --n-clusters (4) crashed inside cluster_molecules_ica. That case now
// falls back to hash initialization safely and is covered by
// tests/test_bugfix_ica_fallback.cpp; these tables keep 4 genes so their
// gene-dependent output expectations stay unchanged.
// with_z adds an exactly per-clump constant z (0/4/8/12) so the table is 3D.
// n_unassigned_per_clump labels that many molecules per clump as "0".
// wide=true places the clumps in a single row along x (~600 vs ~4 in y): the
// HTML reports rasterize 6000px-wide PNGs whose height follows the data aspect
// ratio, so a wide layout keeps those renders fast under -O0.
// collinear_last_clump pins the last clump to a single horizontal line
// (y exactly constant) so at least one resulting cell is degenerate and the
// area/density/elongation fallbacks in the cell-stats code get exercised.
MoleculeTable make_clumped_table(int per_clump, bool with_z,
                                 int n_unassigned_per_clump, unsigned seed = 42,
                                 bool wide = false,
                                 bool collinear_last_clump = false) {
    MoleculeTable t;
    static const char* kGenes[] = {"GeneA", "GeneB", "GeneC", "GeneD"};
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> jit_x(-2.0, 2.0);
    std::uniform_real_distribution<double> jit_y(wide ? -1.0 : -2.0,
                                                 wide ? 1.0 : 2.0);

    const int n_clumps = 4;
    t.x.reserve(static_cast<size_t>(n_clumps) * per_clump);
    t.y.reserve(static_cast<size_t>(n_clumps) * per_clump);
    t.gene.reserve(static_cast<size_t>(n_clumps) * per_clump);
    t.cell.reserve(static_cast<size_t>(n_clumps) * per_clump);
    if (with_z) t.z.reserve(static_cast<size_t>(n_clumps) * per_clump);

    for (int c = 0; c < n_clumps; ++c) {
        const double cx = wide ? (10.0 + 200.0 * c) : (10.0 + 20.0 * (c % 2));
        const double cy = wide ? 10.0 : (10.0 + 20.0 * (c / 2));
        const bool collinear = collinear_last_clump && c == n_clumps - 1;
        for (int i = 0; i < per_clump; ++i) {
            const int idx = t.size();
            t.x.push_back(cx + jit_x(rng));
            t.y.push_back(collinear ? cy : (cy + jit_y(rng)));
            if (with_z) t.z.push_back(4.0 * c);
            t.gene.push_back(kGenes[idx % 4]);
            t.cell.push_back(i < n_unassigned_per_clump
                                 ? "0"
                                 : ("cell" + std::to_string(c + 1)));
        }
    }
    return t;
}

fs::path write_table_csv(const TempDir& tmp, const std::string& name,
                         const MoleculeTable& t, bool with_cell_col) {
    const fs::path p = tmp.path / name;
    std::ofstream f(p);
    f << "x,y,gene";
    if (with_cell_col) f << ",cell_id";
    if (!t.z.empty()) f << ",z";
    f << "\n";
    for (int i = 0; i < t.size(); ++i) {
        f << t.x[i] << "," << t.y[i] << "," << t.gene[i];
        if (with_cell_col) f << "," << t.cell[i];
        if (!t.z.empty()) f << "," << t.z[i];
        f << "\n";
    }
    return p;
}

fs::path write_2d_csv(const TempDir& tmp, int per_clump = 50,
                      const std::string& name = "mols.csv", bool wide = false,
                      bool collinear_last = false) {
    return write_table_csv(tmp, name,
                           make_clumped_table(per_clump, false, 5, 42, wide,
                                              collinear_last),
                           /*with_cell_col=*/true);
}

// Wide-aspect table for tests whose pipeline renders 6000px HTML PNGs.
fs::path write_wide_csv(const TempDir& tmp, int per_clump = 40,
                        const std::string& name = "mols_wide.csv") {
    return write_2d_csv(tmp, per_clump, name, /*wide=*/true);
}

// Binary TIFF mask whose foreground pixels cover every molecule pixel
// (dilated by 1), so each clump becomes one labelled component.
fs::path write_mask_tiff(const fs::path& path, const MoleculeTable& t,
                         uint32_t size = 64) {
    std::vector<uint8_t> px(static_cast<size_t>(size) * size, 0);
    for (int i = 0; i < t.size(); ++i) {
        const int col = static_cast<int>(std::round(t.x[i])) - 1;
        const int row = static_cast<int>(std::round(t.y[i])) - 1;
        for (int dr = -1; dr <= 1; ++dr) {
            for (int dc = -1; dc <= 1; ++dc) {
                const int rr = row + dr;
                const int cc = col + dc;
                if (rr >= 0 && cc >= 0 && rr < static_cast<int>(size) &&
                    cc < static_cast<int>(size)) {
                    px[static_cast<size_t>(rr) * size + cc] = 1;
                }
            }
        }
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

// Boundary polygon that overlaps the molecule bounding box but contains no
// molecules (the clumps occupy [8,12] and [28,32] on both axes).
fs::path write_empty_boundary_csv(const TempDir& tmp,
                                  const std::string& name = "bounds.csv") {
    const fs::path p = tmp.path / name;
    std::ofstream f(p);
    f << "cell_id,vertex_x,vertex_y\n"
      << "cell1,16,16\n"
      << "cell1,16,24\n"
      << "cell1,24,24\n"
      << "cell1,24,16\n";
    return p;
}

fs::path write_text(const TempDir& tmp, const std::string& name,
                    const std::string& content) {
    const fs::path p = tmp.path / name;
    std::ofstream f(p);
    f << content;
    return p;
}

void write_transcripts_parquet(const fs::path& path, const MoleculeTable& t) {
    arrow::DoubleBuilder xb;
    arrow::DoubleBuilder yb;
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

    auto schema = arrow::schema({
        arrow::field("x", arrow::float64()),
        arrow::field("y", arrow::float64()),
        arrow::field("gene", arrow::utf8()),
    });
    auto table = arrow::Table::Make(schema, {xa, ya, ga});

    auto sink_res = arrow::io::FileOutputStream::Open(path.string());
    ASSERT_TRUE(sink_res.ok()) << sink_res.status().ToString();
    auto sink = *sink_res;
    auto st = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), sink, 3);
    EXPECT_TRUE(st.ok()) << st.ToString();
    EXPECT_TRUE(sink->Close().ok());
}

// A Xenium experiment directory: manifest + transcripts.parquet.
fs::path make_xenium_dir(const TempDir& tmp, const MoleculeTable& t,
                         const std::string& dir_name = "xen") {
    const fs::path d = tmp.path / dir_name;
    fs::create_directories(d);
    write_text(tmp, dir_name + "/experiment.xenium",
               "run_folder,panel\nsynthetic,none\n");
    write_transcripts_parquet(d / "transcripts.parquet", t);
    return d / "experiment.xenium";
}

// ---------------------------------------------------------------------------
// Output-file assertions
// ---------------------------------------------------------------------------

bool file_nonempty(const fs::path& p) {
    std::error_code ec;
    return fs::is_regular_file(p, ec) && fs::file_size(p, ec) > 0;
}

void expect_files(const fs::path& dir, const std::vector<std::string>& names) {
    for (const auto& name : names) {
        EXPECT_TRUE(file_nonempty(dir / name))
            << "expected non-empty output file: " << (dir / name);
    }
}

std::string first_line(const fs::path& p) {
    std::ifstream f(p);
    std::string line;
    std::getline(f, line);
    return line;
}

int count_lines(const fs::path& p) {
    std::ifstream f(p);
    int n = 0;
    std::string line;
    while (std::getline(f, line)) ++n;
    return n;
}

// Extract "key = value" from a TOML dump and return the value as string.
std::string toml_value(const std::string& toml, const std::string& key) {
    const std::string prefix = key + " = ";
    size_t pos = toml.find(prefix);
    if (pos == std::string::npos) return "";
    pos += prefix.size();
    size_t end = toml.find('\n', pos);
    return toml.substr(pos, end == std::string::npos ? std::string::npos : end - pos);
}

}  // namespace

// ============================================================================
// --help / --version / CLI11 parse errors
// ============================================================================

TEST(Cov5CliHelp, RootHelpExitsZeroWithSubcommands) {
    TempDir tmp("help_root");
    auto r = run_cli(tmp, "--help");
    EXPECT_EQ(r.exit_code, 0);
    EXPECT_NE(r.out.find("Baysor"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("Run cell segmentation"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("preview"), std::string::npos);
    EXPECT_NE(r.out.find("segfree"), std::string::npos);
}

TEST(Cov5CliHelp, RunSubcommandHelpListsRunOptions) {
    TempDir tmp("help_run");
    auto r = run_cli(tmp, "run --help");
    EXPECT_EQ(r.exit_code, 0);
    EXPECT_NE(r.out.find("--output-style"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("--count-matrix-format"), std::string::npos);
    EXPECT_NE(r.out.find("--nuclei-genes"), std::string::npos);
}

TEST(Cov5CliHelp, VersionFlagPrintsProjectVersion) {
    TempDir tmp("version");
    auto r = run_cli(tmp, "--version");
    EXPECT_EQ(r.exit_code, 0) << r.err;
    EXPECT_EQ(r.out, std::string("baysor ") + BAYSOR_VERSION + "\n");
}

TEST(Cov5CliHelp, SubcommandVersionFlagsPrintBareParseableVersion) {
    TempDir tmp("subcommand_version");
    for (const auto* command : {"run", "preview", "segfree"}) {
        const auto r = run_cli(tmp, std::string(command) + " --version");
        EXPECT_EQ(r.exit_code, 0) << command << ": " << r.err;
        // Sopa parses stdout directly with packaging.version.Version, which
        // accepts the bare version but not the top-level "baysor " prefix.
        EXPECT_EQ(r.out, std::string(BAYSOR_VERSION) + "\n") << command;
    }
}

TEST(Cov5CliParse, RequiresExactlyOneSubcommand) {
    TempDir tmp("parse_no_sub");
    auto r = run_cli(tmp, "");
    EXPECT_NE(r.exit_code, 0);
    EXPECT_FALSE(r.err.empty()) << "expected a CLI11 error on stderr";
}

TEST(Cov5CliParse, UnknownOptionIsRejected) {
    TempDir tmp("parse_unknown");
    auto csv = write_2d_csv(tmp, 20);
    auto r = run_cli(tmp, "run '" + csv.string() + "' -m 10 -s 2.5 --frobnicate");
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.err.find("--frobnicate"), std::string::npos) << r.err;
}

TEST(Cov5CliParse, MissingRequiredCoordinatesIsRejected) {
    TempDir tmp("parse_missing_pos");
    auto r = run_cli(tmp, "run -m 10 -s 2.5");
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.err.find("coordinates"), std::string::npos) << r.err;
}

TEST(Cov5CliParse, NonNumericIntegerOptionIsRejected) {
    TempDir tmp("parse_bad_int");
    auto csv = write_2d_csv(tmp, 20);
    auto r = run_cli(tmp, "run '" + csv.string() + "' -m 10 -s 2.5 --iters abc");
    EXPECT_NE(r.exit_code, 0);
    EXPECT_NE(r.err.find("iters"), std::string::npos) << r.err;
}

// ============================================================================
// Option validation errors from main()'s dispatch
// ============================================================================

TEST(Cov5CliValidate, UnknownOutputStyleFails) {
    TempDir tmp("val_style");
    auto csv = write_2d_csv(tmp, 20);
    auto r = run_cli(tmp, "run '" + csv.string() +
                              "' -m 10 -s 2.5 --output-style xml -o '" +
                              (tmp.path / "seg").string() + "'");
    EXPECT_EQ(r.exit_code, 1);
    EXPECT_NE(r.out.find("Unknown output style: xml"), std::string::npos) << r.out;
}

TEST(Cov5CliValidate, UnknownPolygonFormatFails) {
    TempDir tmp("val_polygon");
    auto csv = write_2d_csv(tmp, 20);
    auto r = run_cli(tmp, "run '" + csv.string() +
                              "' -m 10 -s 2.5 --polygon-format geo -o '" +
                              (tmp.path / "seg").string() + "'");
    EXPECT_NE(r.exit_code, 0);
    const std::string combined = r.out + r.err;
    EXPECT_NE(combined.find("polygon-format"), std::string::npos) << combined;
    EXPECT_FALSE(fs::exists(tmp.path / "seg" / "segmentation_polygons_2d.json"));
}

TEST(Cov5CliValidate, UnknownClusterMethodFails) {
    TempDir tmp("val_cluster");
    auto csv = write_2d_csv(tmp, 20);
    auto r = run_cli(tmp, "run '" + csv.string() +
                              "' -m 10 -s 2.5 --cluster-method kmeans -o '" +
                              (tmp.path / "seg").string() + "'");
    EXPECT_EQ(r.exit_code, 1);
    EXPECT_NE(r.out.find("cluster_method must be one of"), std::string::npos) << r.out;
}

TEST(Cov5CliValidate, MissingBothPriorAndScaleFails) {
    TempDir tmp("val_scale");
    auto csv = write_2d_csv(tmp, 20);
    auto r = run_cli(tmp, "run '" + csv.string() + "' -m 10 -o '" +
                              (tmp.path / "seg").string() + "'");
    EXPECT_EQ(r.exit_code, 1);
    EXPECT_NE(r.out.find("Either prior_segmentation or --scale must be provided."),
              std::string::npos) << r.out;
}

TEST(Cov5CliValidate, MissingConfigFileFails) {
    TempDir tmp("val_config");
    auto csv = write_2d_csv(tmp, 20);
    auto r = run_cli(tmp, "run '" + csv.string() + "' -c '" +
                              (tmp.path / "missing.toml").string() +
                              "' -m 10 -s 2.5");
    EXPECT_EQ(r.exit_code, 1);
    EXPECT_NE(r.out.find("Failed to load config"), std::string::npos) << r.out;
}

// ============================================================================
// Input-file errors (thrown inside cmd_* and caught by main's dispatch try)
// ============================================================================

TEST(Cov5CliInput, MissingCoordinatesFileFails) {
    TempDir tmp("in_missing");
    auto r = run_cli(tmp, "run '" + (tmp.path / "nope.csv").string() +
                              "' -m 10 -s 2.5");
    EXPECT_EQ(r.exit_code, 1);
    // Exact error from load_molecules (arrow_unwrap of the reader open):
    EXPECT_NE(r.out.find("Arrow error: IOError: Failed to open local file"),
              std::string::npos) << r.out;
    EXPECT_NE(r.out.find("No such file or directory"), std::string::npos) << r.out;
}

TEST(Cov5CliInput, MissingPriorMaskFileFails) {
    TempDir tmp("in_mask");
    auto csv = write_2d_csv(tmp, 20);
    auto r = run_cli(tmp, "run '" + csv.string() + "' '" +
                              (tmp.path / "missing_mask.tiff").string() +
                              "' -m 10 -s 2.5");
    EXPECT_EQ(r.exit_code, 1);
    EXPECT_NE(r.out.find("Cannot open TIFF mask file"), std::string::npos) << r.out;
}

TEST(Cov5CliXenium, ManifestWithoutTranscriptsFails) {
    TempDir tmp("xen_missing");
    const fs::path manifest = tmp.path / "experiment.xenium";
    write_text(tmp, "experiment.xenium", "run_folder,panel\nsynthetic,none\n");
    auto r = run_cli(tmp, "run '" + manifest.string() + "' -m 10 -s 2.5");
    EXPECT_EQ(r.exit_code, 1);
    EXPECT_NE(r.out.find("Could not locate transcripts"), std::string::npos) << r.out;
}

// ============================================================================
// Prior-segmentation edge cases
// ============================================================================

TEST(Cov5CliPrior, BoundaryPriorWithoutUsableScaleFails) {
    TempDir tmp("prior_boundary");
    auto csv = write_2d_csv(tmp, 20);
    auto bounds = write_empty_boundary_csv(tmp);
    auto r = run_cli(tmp, "run '" + csv.string() + "' '" + bounds.string() +
                              "' -m 10 -o '" + (tmp.path / "seg").string() + "'");
    EXPECT_EQ(r.exit_code, 1);
    // The boundary contains no molecules, so scale estimation warns and the
    // run aborts before segmentation.
    EXPECT_NE(r.out.find("Could not estimate scale from prior"),
              std::string::npos) << r.out;
    EXPECT_NE(r.out.find("Scale could not be determined"),
              std::string::npos) << r.out;
}

// ============================================================================
// Guarded / failing paths inside cmd_run
// ============================================================================

TEST(Cov5CliRun, NucleiGenesGuardFails) {
    TempDir tmp("run_nuclei");
    auto csv = write_2d_csv(tmp, 20);
    auto r = run_cli(tmp, "run '" + csv.string() +
                              "' -m 10 -s 2.5 --nuclei-genes GeneA -o '" +
                              (tmp.path / "seg").string() + "'");
    EXPECT_EQ(r.exit_code, 1);
    EXPECT_NE(r.out.find("not yet implemented"), std::string::npos) << r.out;
}

TEST(Cov5CliRun, OutputDirectoryCreationFailureIsReported) {
    TempDir tmp("run_mkdir");
    auto csv = write_2d_csv(tmp, 20);
    const fs::path blocker = write_text(tmp, "blocked", "i am a file\n");
    auto r = run_cli(tmp, "run '" + csv.string() + "' -m 10 -s 2.5 -o '" +
                              (blocker / "seg").string() + "'");
    EXPECT_EQ(r.exit_code, 1);
    EXPECT_NE(r.out.find("Could not create output directory"),
              std::string::npos) << r.out;
}

// ============================================================================
// Successful `run` invocations covering the output-style / prior / cluster /
// dimensionality matrix
// ============================================================================

TEST(Cov5CliRun, PriorColumnLegacy2DBundle) {
    TempDir tmp("run_prior_col");
    // Last clump is a perfectly horizontal line: some cell's convex hull is
    // degenerate (area 0), exercising the NaN density/elongation fallbacks.
    auto csv = write_2d_csv(tmp, 50, "mols.csv", /*wide=*/false,
                            /*collinear_last=*/true);
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run '" + csv.string() + "' ':cell_id' -m 10 --iters 12 -o '" +
                              out.string() + "'");
    ASSERT_TRUE(cli_success(r));

    // Prior-aware n_cells_init bookkeeping is logged.
    EXPECT_NE(r.out.find("prior-aware n_cells_init="), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("Segmentation complete"), std::string::npos);

    expect_files(out, {"segmentation.csv", "segmentation_cell_stats.csv",
                       "segmentation_counts.loom", "segmentation_polygons_2d.json",
                       "segmentation_params.dump.toml", "segmentation_log.log"});

    // Per-molecule table: one row per molecule plus header, with the columns
    // that require clustering, NCV colors and assignment confidence.
    EXPECT_EQ(count_lines(out / "segmentation.csv"), 200 + 1);
    const std::string seg_header = first_line(out / "segmentation.csv");
    EXPECT_NE(seg_header.find("cell"), std::string::npos) << seg_header;
    EXPECT_NE(seg_header.find("ncv_color"), std::string::npos) << seg_header;
    EXPECT_NE(seg_header.find("assignment_confidence"), std::string::npos) << seg_header;
    EXPECT_NE(seg_header.find("cluster"), std::string::npos) << seg_header;
    EXPECT_EQ(seg_header.find(",z"), std::string::npos) << seg_header;

    // Cell stats carry the clustered/lifespan/assignment-confidence columns.
    const std::string stats_header = first_line(out / "segmentation_cell_stats.csv");
    EXPECT_NE(stats_header.find("cluster"), std::string::npos) << stats_header;
    EXPECT_NE(stats_header.find("max_cluster_frac"), std::string::npos) << stats_header;
    EXPECT_NE(stats_header.find("lifespan"), std::string::npos) << stats_header;
    EXPECT_NE(stats_header.find("avg_assignment_confidence"), std::string::npos)
        << stats_header;
    EXPECT_GE(count_lines(out / "segmentation_cell_stats.csv"), 3);
    // The collinear clump yields a zero-area cell: density/elongation fall
    // back to NaN. Check the actual columns rather than any "nan" substring:
    // locate density and elongation in the header and require "nan" in both
    // columns of at least one data row.
    {
        std::ifstream stats_f(out / "segmentation_cell_stats.csv");
        std::string header;
        ASSERT_TRUE(std::getline(stats_f, header));
        auto col_index = [&](const std::string& name) -> int {
            std::stringstream hs(header);
            std::string col;
            int idx = 0;
            while (std::getline(hs, col, ',')) {
                if (col == name) return idx;
                ++idx;
            }
            return -1;
        };
        const int density_col = col_index("density");
        const int elongation_col = col_index("elongation");
        ASSERT_GE(density_col, 0) << header;
        ASSERT_GE(elongation_col, 0) << header;
        bool density_nan = false;
        bool elongation_nan = false;
        std::string row;
        while (std::getline(stats_f, row)) {
            std::stringstream rs(row);
            std::string cell;
            std::vector<std::string> cols;
            while (std::getline(rs, cell, ',')) cols.push_back(cell);
            ASSERT_GT(cols.size(),
                      static_cast<size_t>(std::max(density_col, elongation_col)));
            if (cols[density_col] == "nan") density_nan = true;
            if (cols[elongation_col] == "nan") elongation_nan = true;
        }
        EXPECT_TRUE(density_nan) << "no nan density in:\n"
                                 << read_text_file(out / "segmentation_cell_stats.csv");
        EXPECT_TRUE(elongation_nan) << "no nan elongation in:\n"
                                    << read_text_file(out / "segmentation_cell_stats.csv");
    }

    EXPECT_NE(read_text_file(out / "segmentation_polygons_2d.json")
                  .find("FeatureCollection"),
              std::string::npos);

    // Params dump records the prior and a positive scale estimated from it.
    const std::string dump = read_text_file(out / "segmentation_params.dump.toml");
    EXPECT_NE(toml_value(dump, "type"), "") << dump;
    EXPECT_EQ(toml_value(dump, "type"), "\"column\"") << dump;
    EXPECT_GT(std::stod(toml_value(dump, "scale")), 0.0) << dump;
    EXPECT_EQ(toml_value(dump, "iters"), "12") << dump;

    // Dual logger: messages also land in the log file.
    EXPECT_NE(read_text_file(out / "segmentation_log.log")
                  .find("Segmentation complete"),
              std::string::npos);
}

TEST(Cov5CliRun, PriorColumnLegacyGeometryCollectionLegacyHasIntegerCellIds) {
    TempDir tmp("run_legacy_poly");
    // Last clump is a horizontal line: its free-form polygon estimation fails
    // and the fallback keeps it in the polygons file (kharchenkolab/Baysor#165).
    auto csv = write_2d_csv(tmp, 50, "mols_legacy.csv", /*wide=*/false,
                            /*collinear_last=*/true);
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run '" + csv.string() +
                              "' ':cell_id' -m 10 --iters 12 -s 2.5"
                              " --polygon-format GeometryCollectionLegacy -o '" +
                              out.string() + "'");
    ASSERT_TRUE(cli_success(r));

    const std::string json_text = read_text_file(out / "segmentation_polygons_2d.json");
    const auto doc = nlohmann::json::parse(json_text);
    EXPECT_EQ(doc.at("type"), "GeometryCollection");

    // Collect integer polygon ids and the CSV's cell_<n> names; the sets must
    // be identical after stripping the CSV prefix.
    std::set<int> poly_ids;
    for (const auto& geom : doc.at("geometries")) {
        ASSERT_TRUE(geom.at("cell").is_number_integer()) << geom.dump();
        poly_ids.insert(geom.at("cell").get<int>());
    }
    std::set<int> csv_ids;
    {
        std::ifstream f(out / "segmentation.csv");
        std::string header;
        ASSERT_TRUE(std::getline(f, header));
        const auto cols = header.find("cell,");
        ASSERT_NE(cols, std::string::npos) << header;
        const int cell_col = static_cast<int>(std::count(header.begin(),
                                                         header.begin() + cols, ','));
        std::string line;
        while (std::getline(f, line)) {
            std::stringstream ss(line);
            std::string field;
            for (int i = 0; i <= cell_col; ++i) ASSERT_TRUE(std::getline(ss, field, ','));
            if (field.rfind("cell_", 0) == 0) {
                csv_ids.insert(std::stoi(field.substr(5)));
            }
        }
    }
    EXPECT_EQ(poly_ids, csv_ids);
    EXPECT_FALSE(poly_ids.empty());
}

TEST(Cov5CliRun, ParquetStyleWithPlotAndIgnoredFormatWarnings) {
    TempDir tmp("run_parquet_plot");
    auto csv = write_wide_csv(tmp, 40);
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run '" + csv.string() +
                              "' -m 10 -s 2.5 --iters 12 -p"
                              " --output-style parquet"
                              " --cluster-method none"
                              " --polygon-format GeometryCollection"
                              " --count-matrix-format tsv -o '" +
                              out.string() + "'");
    ASSERT_TRUE(cli_success(r));

    // Non-default polygon/count formats are ignored under parquet output.
    EXPECT_NE(r.out.find("--polygon-format is ignored for output style 'parquet'"),
              std::string::npos) << r.out;
    EXPECT_NE(r.out.find("--count-matrix-format is ignored for output style 'parquet'"),
              std::string::npos) << r.out;
    EXPECT_NE(r.out.find("Generating HTML run report"), std::string::npos);

    expect_files(out, {"molecules.parquet", "cells.parquet", "cell_boundaries.parquet",
                       "feature_matrix.h5", "diagnostic_report.html",
                       "segmentation_plot.html", "run_params.toml", "run.log"});
    EXPECT_NE(read_text_file(out / "diagnostic_report.html").find("html"),
              std::string::npos);
    EXPECT_NE(read_text_file(out / "segmentation_plot.html").find("html"),
              std::string::npos);

    // No legacy-style files were written.
    EXPECT_FALSE(fs::exists(out / "segmentation.csv"));
    EXPECT_FALSE(fs::exists(out / "segmentation_counts.loom"));
}

TEST(Cov5CliRun, ConfigFileWithCliOverride3DAndLouvain) {
    TempDir tmp("run_config3d");
    auto table = make_clumped_table(30, /*with_z=*/true, 5);
    auto csv = write_table_csv(tmp, "mols3d.csv", table, /*with_cell_col=*/true);
    const fs::path cfg = write_text(tmp, "baysor.toml",
                                    "[molecules]\n"
                                    "min_molecules_per_cell = 10\n"
                                    "\n"
                                    "[segmentation]\n"
                                    "scale = 6.0\n"
                                    "cluster_method = \"louvain\"\n"
                                    "iters = 400\n");
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run -c '" + cfg.string() + "' '" + csv.string() +
                              "' --iters 12 --polygon-format GeometryCollection -o '" +
                              out.string() + "'");
    ASSERT_TRUE(cli_success(r));

    EXPECT_NE(r.out.find("Loaded config from"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("Louvain on NCV kNN graph"), std::string::npos) << r.out;

    expect_files(out, {"segmentation.csv", "segmentation_cell_stats.csv",
                       "segmentation_counts.loom", "segmentation_polygons_3d.json",
                       "segmentation_params.dump.toml"});

    // 3D input produces a z column in the cell stats.
    const std::string stats_header = first_line(out / "segmentation_cell_stats.csv");
    EXPECT_NE(stats_header.find(",z"), std::string::npos) << stats_header;

    // Config supplies defaults; the explicit --iters flag wins.
    const std::string dump = read_text_file(out / "segmentation_params.dump.toml");
    EXPECT_EQ(toml_value(dump, "cluster_method"), "\"louvain\"") << dump;
    EXPECT_EQ(toml_value(dump, "min_molecules_per_cell"), "10") << dump;
    EXPECT_EQ(toml_value(dump, "iters"), "12") << dump;
    EXPECT_EQ(toml_value(dump, "iters"), "12");
    EXPECT_EQ(dump.find("iters = 400"), std::string::npos) << dump;
    EXPECT_EQ(std::stod(toml_value(dump, "scale")), 6.0) << dump;
}

TEST(Cov5CliRun, Leiden3DParquetBundle) {
    TempDir tmp("run_leiden3d");
    auto table = make_clumped_table(30, /*with_z=*/true, 5);
    auto csv = write_table_csv(tmp, "mols3d.csv", table, /*with_cell_col=*/true);
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run '" + csv.string() +
                              "' -m 10 -s 2.5 --iters 12"
                              " --cluster-method leiden --output-style parquet -o '" +
                              out.string() + "'");
    ASSERT_TRUE(cli_success(r));

    EXPECT_NE(r.out.find("Leiden on NCV kNN graph"), std::string::npos) << r.out;
    expect_files(out, {"molecules.parquet", "cells.parquet", "feature_matrix.h5",
                       "run_params.toml", "run.log"});
    // 3D input writes both GeoParquet boundary files: the combined 2-D
    // projection and the per-layer 3-D stack.
    expect_files(out, {"cell_boundaries.parquet", "cell_boundaries_3d.parquet"});
    const std::string dump = read_text_file(out / "run_params.toml");
    EXPECT_EQ(toml_value(dump, "cluster_method"), "\"leiden\"") << dump;
}

TEST(Cov5CliRun, PlotHtmlWriteFailuresAreReported) {
    TempDir tmp("run_plot_fail");
    auto csv = write_wide_csv(tmp, 20);
    const fs::path out = tmp.path / "seg";
    fs::create_directories(out / "diagnostic_report.html");
    fs::create_directories(out / "segmentation_plot.html");

    auto r = run_cli(tmp, "run '" + csv.string() +
                              "' -m 10 -s 2.5 --iters 10 -p"
                              " --cluster-method none --skip-ncv-color -o '" +
                              out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    EXPECT_NE(r.out.find("Could not write diagnostic report"),
              std::string::npos) << r.out;
    EXPECT_NE(r.out.find("Could not write segmentation plot"),
              std::string::npos) << r.out;
}

TEST(Cov5CliRun, ImagePriorNoClusteringTsvWithoutPolygons) {
    TempDir tmp("run_image");
    auto table = make_clumped_table(50, /*with_z=*/false, 0);
    auto csv = write_table_csv(tmp, "mols.csv", table, /*with_cell_col=*/true);
    auto mask = write_mask_tiff(tmp.path / "mask.tif", table);
    const fs::path out = tmp.path / "seg";

    auto r = run_cli(tmp, "run '" + csv.string() + "' '" + mask.string() +
                              "' -m 10 --iters 12 --cluster-method none"
                              " --count-matrix-format tsv --polygon-format none"
                              " --skip-ncv-color -o '" +
                              out.string() + "'");
    ASSERT_TRUE(cli_success(r));

    expect_files(out, {"segmentation.csv", "segmentation_cell_stats.csv",
                       "segmentation_counts.tsv", "segmentation_params.dump.toml"});

    // tsv count matrix: header row keyed by gene, one row per gene.
    const std::string tsv = read_text_file(out / "segmentation_counts.tsv");
    EXPECT_EQ(tsv.find("gene\tcell_1"), 0u) << tsv.substr(0, 80);
    EXPECT_NE(tsv.find("GeneA"), std::string::npos) << tsv;

    // polygon-format none + legacy style: no polygon file at all.
    EXPECT_FALSE(fs::exists(out / "segmentation_polygons_2d.json"));

    // cluster-method none: no cluster / max_cluster_frac columns in stats.
    const std::string stats_header = first_line(out / "segmentation_cell_stats.csv");
    EXPECT_NE(stats_header.find("n_transcripts"), std::string::npos) << stats_header;
    EXPECT_EQ(stats_header.find("max_cluster_frac"), std::string::npos) << stats_header;
    EXPECT_EQ(stats_header.find(",cluster"), std::string::npos) << stats_header;

    // --skip-ncv-color drops the ncv_color column from the molecule table.
    const std::string seg_header = first_line(out / "segmentation.csv");
    EXPECT_EQ(seg_header.find("ncv_color"), std::string::npos) << seg_header;

    const std::string dump = read_text_file(out / "segmentation_params.dump.toml");
    EXPECT_EQ(toml_value(dump, "type"), "\"image\"") << dump;
    EXPECT_GT(std::stod(toml_value(dump, "scale")), 0.0) << dump;
}

// ============================================================================
// preview
// ============================================================================

TEST(Cov5CliPreview, WritesHtmlReport) {
    TempDir tmp("preview_ok");
    auto csv = write_wide_csv(tmp, 40);
    const fs::path out = tmp.path / "preview.html";

    auto r = run_cli(tmp, "preview '" + csv.string() + "' -m 10 -o '" +
                              out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    EXPECT_NE(r.out.find("Preview saved to"), std::string::npos) << r.out;
    ASSERT_TRUE(file_nonempty(out));
    const std::string html = read_text_file(out);
    EXPECT_NE(html.find("html"), std::string::npos);
}

TEST(Cov5CliPreview, UnwritableOutputFails) {
    TempDir tmp("preview_fail");
    auto csv = write_wide_csv(tmp, 20);
    const fs::path out = tmp.path / "no_such_dir" / "preview.html";

    auto r = run_cli(tmp, "preview '" + csv.string() + "' -m 10 -o '" +
                              out.string() + "'");
    EXPECT_EQ(r.exit_code, 1);
    EXPECT_NE(r.out.find("Could not write to"), std::string::npos) << r.out;
}

// ============================================================================
// segfree
// ============================================================================

TEST(Cov5CliSegfree, InfersDefaultKAndWritesLoom) {
    TempDir tmp("segfree_default");
    auto csv = write_2d_csv(tmp, 50);
    const fs::path out = tmp.path / "ncvs.loom";

    auto r = run_cli(tmp, "segfree '" + csv.string() + "' -m 10 -o '" +
                              out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    // default_param_value("composition_neighborhood", min_molecules_per_cell=10,
    // n_genes=4) == 10
    EXPECT_NE(r.out.find("Using k=10 neighbors for NCV composition"),
              std::string::npos) << r.out;
    EXPECT_TRUE(file_nonempty(out)) << out;
}

TEST(Cov5CliSegfree, ExplicitKWithXeniumManifestInput) {
    TempDir tmp("segfree_xenium");
    auto table = make_clumped_table(50, /*with_z=*/false, 0);
    const fs::path manifest = make_xenium_dir(tmp, table);
    const fs::path out = tmp.path / "ncvs.loom";

    auto r = run_cli(tmp, "segfree '" + manifest.string() + "' -m 10 -k 12 -o '" +
                              out.string() + "'");
    ASSERT_TRUE(cli_success(r));
    EXPECT_NE(r.out.find("Using k=12 neighbors for NCV composition"),
              std::string::npos) << r.out;
    EXPECT_TRUE(file_nonempty(out)) << out;
}

#endif  // BAYSOR_CLI_PATH && !defined(_WIN32)
