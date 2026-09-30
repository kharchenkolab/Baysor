// Determinism regression tests for the graph/segmentation pipeline.
//
// Root cause covered here (DET-RUNID): CGAL's finite_edges iterator emits
// each edge from whichever adjacent face has the lower *heap address*
// (CGAL Triangulation_ds_iterators_2.h, associated_edge() compares raw
// Face_handle pointers). The edge emission order therefore depends on the
// heap layout, which shifts with the length of the `-o` output path and with
// allocator history (e.g. concurrent parquet decoding). That order flows into
// the CSR adjacency lists and into floating-point summation order of the
// stochastic MRF E-step, so 1-thread runs used to depend on the output path
// length.
//
// The fix canonicalises the edge stream (sorted by (src, dst)) before dedup
// and CSR construction, making the graph a pure function of the edge set.
//
// Tests:
//   1. AdjacencyEdgesAreCanonicallySorted — direct contract on the fix.
//   2. AdjacencyOrderStableUnderHeapPerturbation — the behavioural guard:
//      the same points must give byte-identical adjacency lists regardless of
//      preceding heap activity.
//   3. SameAssignmentsForDifferentOutputPathLengthsAtOneThread (POSIX/CLI) —
//      end-to-end: the same tiny dataset segmented at OMP_NUM_THREADS=1 with
//      two output paths of different lengths must produce identical
//      assignments.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "baysor/processing/data_processing/triangulation.h"

namespace {

// Deterministic pseudo-random 2D point cloud (no global RNG, no <random>
// machinery that could interact with other tests).
Eigen::MatrixXd make_det_points(int n) {
    Eigen::MatrixXd pts(2, n);
    std::uint64_t s = 12345;
    auto next = [&s]() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>((s >> 16) & 0xFFFFFF) / 16777216.0;
    };
    for (int i = 0; i < n; ++i) {
        pts(0, i) = next() * 100.0;
        pts(1, i) = next() * 100.0;
    }
    return pts;
}

bool same_edges(const baysor::AdjacencyResult& a, const baysor::AdjacencyResult& b) {
    if (a.edge_src != b.edge_src) return false;
    if (a.edge_dst != b.edge_dst) return false;
    if (a.edge_dists.size() != b.edge_dists.size()) return false;
    for (size_t i = 0; i < a.edge_dists.size(); ++i) {
        // Bitwise equality: the order must not even perturb the last bits.
        std::uint64_t ba, bb;
        std::memcpy(&ba, &a.edge_dists[i], sizeof(ba));
        std::memcpy(&bb, &b.edge_dists[i], sizeof(bb));
        if (ba != bb) return false;
    }
    return true;
}

}  // namespace

TEST(DetGraph, AdjacencyEdgesAreCanonicallySorted) {
    const Eigen::MatrixXd pts = make_det_points(600);
    auto res = baysor::adjacency_list(pts, /*filter=*/false);
    ASSERT_GT(res.edge_src.size(), 0u);

    for (size_t i = 1; i < res.edge_src.size(); ++i) {
        const auto prev = std::minmax(res.edge_src[i - 1], res.edge_dst[i - 1]);
        const auto cur = std::minmax(res.edge_src[i], res.edge_dst[i]);
        ASSERT_TRUE(prev.first < cur.first ||
                    (prev.first == cur.first && prev.second < cur.second))
            << "edges not canonically sorted at index " << i;
    }
}

TEST(DetGraph, AdjacencyOrderStableUnderHeapPerturbation) {
    const Eigen::MatrixXd pts = make_det_points(600);
    const auto ref = baysor::adjacency_list(pts, /*filter=*/false);
    ASSERT_GT(ref.edge_src.size(), 0u);

    // Allocate and hold differently sized/fragmented chunks between builds so
    // subsequent allocations (the triangulation's internal blocks) land in
    // different places relative to each other — the exact condition under
    // which the unsorted, address-dependent edge emission order changed.
    for (int perturb = 0; perturb < 6; ++perturb) {
        std::vector<std::unique_ptr<char[]>> held;
        for (int k = 0; k <= perturb; ++k) {
            // Sizes spanning the CGAL block allocations (blocks grow from
            // ~1 KiB upwards) plus a couple of larger ones.
            const size_t sz = 1024 + static_cast<size_t>(k) * 1536 + 64;
            held.push_back(std::make_unique<char[]>(sz));
            std::memset(held.back().get(), 0x5A, sz);
        }
        auto got = baysor::adjacency_list(pts, /*filter=*/false);
        EXPECT_TRUE(same_edges(ref, got))
            << "adjacency list changed under heap perturbation " << perturb;
        // held freed at end of iteration
    }
}

// ---------------------------------------------------------------------------
// End-to-end CLI check (POSIX only): two output path lengths, 1 thread.
// ---------------------------------------------------------------------------

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)

#include <sys/wait.h>
#include <unistd.h>

namespace {

namespace fs = std::filesystem;

// Minimal clumped dataset writer (independent of test_cov_cli.cpp helpers).
fs::path write_tiny_csv(const fs::path& dir) {
    const fs::path p = dir / "mols.csv";
    std::ofstream f(p);
    f << "x,y,gene,cell_id\n";
    const char* genes[] = {"GeneA", "GeneB", "GeneC", "GeneD"};
    std::uint64_t s = 42;
    auto next = [&s]() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>((s >> 16) & 0xFFFFFF) / 16777216.0;
    };
    int idx = 0;
    for (int c = 0; c < 4; ++c) {
        const double cx = 10.0 + 20.0 * (c % 2);
        const double cy = 10.0 + 20.0 * (c / 2);
        for (int i = 0; i < 100; ++i, ++idx) {
            const double x = cx + (next() * 4.0 - 2.0);
            const double y = cy + (next() * 4.0 - 2.0);
            f << x << "," << y << "," << genes[idx % 4] << ","
              << (i >= 5 ? ("cell" + std::to_string(c + 1)) : "0") << "\n";
        }
    }
    return p;
}

int run_cli_raw(const std::string& cmd) {
    const int status = std::system(cmd.c_str());
    if (status >= 0 && WIFEXITED(status)) return WEXITSTATUS(status);
    return -1;
}

std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

}  // namespace

TEST(DetCli, SameAssignmentsForDifferentOutputPathLengthsAtOneThread) {
    const fs::path root =
        fs::temp_directory_path() / "baysor_det_runid_pathlen";
    std::error_code ec;
    fs::remove_all(root, ec);
    fs::create_directories(root);

    const fs::path csv = write_tiny_csv(root);

    // Same command, same data, same environment (OMP_NUM_THREADS=1); only
    // the length of the -o path differs (by 16 characters).
    const fs::path out_short = root / "o1" / "seg";
    const fs::path out_long = root / "o1_padding_16_chars_here" / "seg";

    const std::string common =
        "'" + std::string(BAYSOR_CLI_PATH) + "' run '" + csv.string() +
        "' -m 10 -s 2.5 --iters 30 --cluster-method none -o ";

    const int rc1 = run_cli_raw("env OMP_NUM_THREADS=1 " + common + "'" +
                                out_short.string() + "' > '" +
                                (root / "r1.out").string() + "' 2> '" +
                                (root / "r1.err").string() + "'");
    const int rc2 = run_cli_raw("env OMP_NUM_THREADS=1 " + common + "'" +
                                out_long.string() + "' > '" +
                                (root / "r2.out").string() + "' 2> '" +
                                (root / "r2.err").string() + "'");
    ASSERT_EQ(rc1, 0) << read_file(root / "r1.out") << read_file(root / "r1.err");
    ASSERT_EQ(rc2, 0) << read_file(root / "r2.out") << read_file(root / "r2.err");

    const fs::path seg1 = out_short / "segmentation.csv";
    const fs::path seg2 = out_long / "segmentation.csv";
    ASSERT_TRUE(fs::exists(seg1)) << seg1;
    ASSERT_TRUE(fs::exists(seg2)) << seg2;

    const std::string s1 = read_file(seg1);
    const std::string s2 = read_file(seg2);
    ASSERT_GT(s1.size(), 0u);
    EXPECT_EQ(s1, s2)
        << "segmentation depends on the output path length at 1 thread";

    fs::remove_all(root, ec);
}

#else

TEST(DetCli, SameAssignmentsForDifferentOutputPathLengthsAtOneThread) {
    GTEST_SKIP() << "CLI subprocess tests require POSIX and BAYSOR_CLI_PATH";
}

#endif
