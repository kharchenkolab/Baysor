// Regression test for kharchenkolab/Baysor#117:
// `--prior-segmentation-confidence 1` must keep every prior cell together.
//
// The invariant is: all molecules carrying the same prior-segmentation label
// end up in one final cell (cell ids may be renamed and cells may grow, but a
// prior cell must never be split across several final cells or partly
// dropped to noise). The reporter saw it violated on CosMx data, where
// duplicated transcripts at identical coordinates are common; overlapping
// prior cells are enough to break it even without duplicates.
//
// The fixture is two overlapping prior cells (a 10x4 grid at (0,0) and the
// same grid shifted by 1.0 in x, `--scale 1.2`), optionally with exact
// duplicate coordinates. Before the fix the E-step assigns a boundary
// molecule to whichever adjacent component locally dominates, so both prior
// cells end up spread over both final cells. The assertion below is the
// invariant itself, so it stays valid for any future algorithm change.
//
// Suite name is prefixed with Bug117 to avoid clashes with the other files.

#include <gtest/gtest.h>

#include "test_cov_helpers.h"

#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

using baysor_test::TempDir;
using baysor_test::cli::CliResult;
using baysor_test::cli::read_text_file;
using baysor_test::cli::run_cli;

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)

// Two overlapping prior cells; `dups` adds exact duplicates of two molecules
// per cell (same coordinates and label, a different gene), matching the CosMx
// duplicate-transcript pattern from the issue.
std::string make_overlapping_prior_csv(const TempDir& tmp,
                                       const std::string& name, bool dups) {
    const char* genes[] = {"GeneA", "GeneB", "GeneC", "GeneD"};
    std::ostringstream f;
    f << "x,y,gene,cell_id\n";
    int k = 0;
    std::vector<std::string> dup_rows;
    for (int cell = 1; cell <= 2; ++cell) {
        const double cx = (cell == 1) ? 0.0 : 1.0;
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 4; ++j) {
                const double x = cx - 0.9 + 0.2 * i;
                const double y = -0.9 + 0.6 * j;
                f << x << ',' << y << ',' << genes[k % 4] << ",cell" << cell << '\n';
                if (dups && (i == 3 || i == 7) && j == 1) {
                    std::ostringstream d;
                    d << x << ',' << y << ',' << genes[(k + 1) % 4] << ",cell"
                      << cell << '\n';
                    dup_rows.push_back(d.str());
                }
                ++k;
            }
        }
    }
    for (const auto& row : dup_rows) f << row;

    const std::string p = tmp.file(name);
    std::ofstream out(p);
    out << f.str();
    return p;
}

// Prior label per molecule in input order (the FASTA-like fixture above).
std::vector<std::string> prior_labels(bool dups) {
    std::vector<std::string> labels;
    for (int cell = 1; cell <= 2; ++cell) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 4; ++j) {
                labels.push_back("cell" + std::to_string(cell));
            }
        }
    }
    if (dups) labels.push_back("cell1"), labels.push_back("cell1");
    if (dups) labels.push_back("cell2"), labels.push_back("cell2");
    return labels;
}

// Reads the "cell" column of a legacy segmentation.csv, in row order.
std::vector<std::string> read_cell_column(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    EXPECT_TRUE(static_cast<bool>(std::getline(f, line))) << path;
    const std::size_t comma = line.find(',');
    EXPECT_EQ(line.substr(0, comma), "cell");
    std::vector<std::string> cells;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        cells.push_back(line.substr(0, line.find(',')));
    }
    return cells;
}

void expect_invariant_holds(const TempDir& tmp, const std::string& csv,
                            bool dups) {
    const std::vector<std::string> labels = prior_labels(dups);

    // The invariant must hold at any thread count, and the result must not
    // depend on it (chunk-keyed E-step RNG streams plus the deterministic
    // projection below).
    std::vector<std::vector<std::string>> cells_per_threads;
    for (const int threads : {1, 2}) {
        const std::string out =
            (tmp.path / ("seg_t" + std::to_string(threads))).string();
        CliResult r = run_cli(tmp, "run '" + csv +
                                       "' ':cell_id' -s 1.2 -m 5 -t " +
                                       std::to_string(threads) +
                                       " --iters 12 --cluster-method none"
                                       " --skip-ncv-color"
                                       " --prior-segmentation-confidence 1"
                                       " -o '" + out + "'");
        ASSERT_EQ(r.exit_code, 0) << r.out << "\n" << r.err;
        cells_per_threads.push_back(
            read_cell_column(out + "/segmentation.csv"));
        ASSERT_EQ(cells_per_threads.back().size(), labels.size());
    }
    EXPECT_EQ(cells_per_threads[0], cells_per_threads[1])
        << "segmentation depends on the thread count at confidence 1";
    const std::vector<std::string>& cells = cells_per_threads[0];

    std::map<std::string, std::set<std::string>> cells_per_label;
    std::map<std::string, int> noise_per_label;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        if (cells[i] == "0") {
            ++noise_per_label[labels[i]];
        } else {
            cells_per_label[labels[i]].insert(cells[i]);
        }
    }

    for (const auto& [label, c] : cells_per_label) {
        EXPECT_EQ(c.size(), 1u) << "prior cell " << label
                                << " was split across " << c.size()
                                << " final cells";
    }
    for (const auto& [label, n] : noise_per_label) {
        EXPECT_EQ(n, 0) << "prior cell " << label << " has " << n
                        << " molecules in noise";
    }
    // Every prior label must appear in the maps (no label silently missing).
    EXPECT_EQ(cells_per_label.size(), static_cast<std::size_t>(2));
}

#endif  // !defined(_WIN32) && defined(BAYSOR_CLI_PATH)

}  // namespace

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)

TEST(Bug117PriorConsistency, ConfidenceOneKeepsOverlappingPriorCellsTogether) {
    TempDir tmp("bug117_prior");
    const std::string csv = make_overlapping_prior_csv(tmp, "mols.csv", /*dups=*/false);
    expect_invariant_holds(tmp, csv, /*dups=*/false);
}

TEST(Bug117PriorConsistency, ConfidenceOneKeepsDuplicatedMoleculesTogether) {
    TempDir tmp("bug117_prior_dup");
    const std::string csv = make_overlapping_prior_csv(tmp, "mols.csv", /*dups=*/true);
    expect_invariant_holds(tmp, csv, /*dups=*/true);
}

#else

TEST(Bug117PriorConsistency, RequiresPosixAndCliPath) {
    GTEST_SKIP() << "requires POSIX and BAYSOR_CLI_PATH";
}

#endif
