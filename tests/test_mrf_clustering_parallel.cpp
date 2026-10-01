// MRF molecule clustering (cluster_molecules_on_mrf) runs its E-step,
// convergence statistics and M-step in parallel. All three must give bitwise
// identical results at any thread count: the E-step is per molecule, the
// convergence statistics are exact max / count reductions over fixed chunks,
// and the M-step accumulates every gene in molecule order.

#include <gtest/gtest.h>

#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/utils/thread_pool.h"

#include <random>
#include <vector>

namespace {

class PoolSizeGuard {
public:
    explicit PoolSizeGuard(int n) : old_(baysor::thread_pool_size()) {
        baysor::set_thread_pool_size(n);
    }
    ~PoolSizeGuard() { baysor::set_thread_pool_size(old_); }

private:
    int old_;
};

struct MrfInput {
    std::vector<int> genes;  // 1-based; 0 = unknown gene
    std::vector<double> confidence;
    baysor::AdjList adj;
};

// Molecules on a side x side grid with 4-neighbour edges. Two spatial domains
// with different gene preferences, a few unknown genes, random confidences.
MrfInput make_grid_input(int side, int n_genes, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const int n = side * side;
    MrfInput in;
    in.genes.resize(n);
    in.confidence.resize(n);
    for (int y = 0; y < side; ++y) {
        for (int x = 0; x < side; ++x) {
            const int i = y * side + x;
            const bool left = x < side / 2;
            int g = 1 + static_cast<int>(unif(rng) * (n_genes / 2));
            if (!left) g += n_genes / 2;
            if (unif(rng) < 0.2) g = 1 + static_cast<int>(unif(rng) * n_genes);
            if (unif(rng) < 0.01) g = 0;
            in.genes[i] = std::min(g, n_genes);
            in.confidence[i] = 0.3 + 0.7 * unif(rng);
        }
    }
    in.genes[0] = n_genes;  // make sure the largest gene id is present
    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int y = 0; y < side; ++y) {
        for (int x = 0; x < side; ++x) {
            const int i = y * side + x;
            if (x + 1 < side) { src.push_back(i); dst.push_back(i + 1); wts.push_back(0.5 + unif(rng)); }
            if (y + 1 < side) { src.push_back(i); dst.push_back(i + side); wts.push_back(0.5 + unif(rng)); }
        }
    }
    in.adj = baysor::AdjList::from_edge_list(
        src.data(), dst.data(), wts.data(), static_cast<int>(src.size()), n);
    return in;
}

void expect_bitwise_equal(const baysor::ClusteringResult& a, const baysor::ClusteringResult& b) {
    ASSERT_EQ(a.assignment, b.assignment);
    ASSERT_EQ(a.diffs, b.diffs);
    ASSERT_EQ(a.change_fracs, b.change_fracs);
    ASSERT_EQ(a.assignment_probs.rows(), b.assignment_probs.rows());
    ASSERT_EQ(a.assignment_probs.cols(), b.assignment_probs.cols());
    for (Eigen::Index i = 0; i < a.assignment_probs.size(); ++i)
        ASSERT_EQ(a.assignment_probs.data()[i], b.assignment_probs.data()[i]) << "prob " << i;
    ASSERT_EQ(a.exprs.rows(), b.exprs.rows());
    ASSERT_EQ(a.exprs.cols(), b.exprs.cols());
    for (Eigen::Index i = 0; i < a.exprs.size(); ++i)
        ASSERT_EQ(a.exprs.data()[i], b.exprs.data()[i]) << "expr " << i;
}

} // namespace

TEST(MrfClusteringParallel, ResultsDoNotDependOnThreadCount) {
    // 80 x 80 = 6,400 molecules: 13 E-step chunks of 512 molecules.
    const MrfInput in = make_grid_input(80, 40, 17);
    baysor::ClusteringResult ref;
    {
        PoolSizeGuard pool(1);
        ref = baysor::cluster_molecules_on_mrf(in.genes, in.adj, in.confidence, 3, 0.01, 1.0, 300, false);
    }
    ASSERT_FALSE(ref.diffs.empty());
    ASSERT_EQ(ref.diffs.size(), ref.change_fracs.size());
    for (int threads : {2, 3, 8}) {
        PoolSizeGuard pool(threads);
        auto res = baysor::cluster_molecules_on_mrf(in.genes, in.adj, in.confidence, 3, 0.01, 1.0, 300, false);
        SCOPED_TRACE(threads);
        expect_bitwise_equal(ref, res);
    }
}

TEST(MrfClusteringParallel, ConvergenceTraceIsConsistent) {
    const MrfInput in = make_grid_input(40, 20, 3);
    PoolSizeGuard pool(4);
    auto res = baysor::cluster_molecules_on_mrf(in.genes, in.adj, in.confidence, 2, 0.01, 1.0, 500, false);
    ASSERT_FALSE(res.diffs.empty());
    for (size_t t = 0; t < res.diffs.size(); ++t) {
        EXPECT_GE(res.diffs[t], 0.0);
        EXPECT_LE(res.diffs[t], 1.0);
        EXPECT_GE(res.change_fracs[t], 0.0);
        EXPECT_LE(res.change_fracs[t], 1.0);
        // A molecule counts as changed only above 1e-7, so no change at all
        // means the maximum is at most 1e-7.
        if (res.change_fracs[t] == 0.0) EXPECT_LE(res.diffs[t], 1e-7);
    }
    // Converged: the last 21 maxima are all below tol.
    if (static_cast<int>(res.diffs.size()) < 500) {
        for (size_t t = res.diffs.size() - 21; t < res.diffs.size(); ++t)
            EXPECT_LT(res.diffs[t], 0.01);
    }
}
