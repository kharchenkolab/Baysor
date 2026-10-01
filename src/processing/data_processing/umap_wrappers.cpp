#include "baysor/processing/data_processing/umap_wrappers.h"
#include "baysor/utils/thread_pool.h"

// Route the FetchContent dependencies' parallel loops through the Baysor
// thread pool instead of OpenMP or per-call std::thread spawns. These macros
// must be defined before any subpar/umappp header is included; the hook
// functions themselves are declared in baysor/utils/thread_pool.h above.
#define SUBPAR_CUSTOM_PARALLELIZE_RANGE baysor::subpar_parallelize_range
#define SUBPAR_CUSTOM_PARALLELIZE_RANGE_NOTHROW baysor::subpar_parallelize_range
#define SUBPAR_CUSTOM_PARALLELIZE_SIMPLE baysor::subpar_parallelize_simple
#define SUBPAR_CUSTOM_PARALLELIZE_SIMPLE_NOTHROW baysor::subpar_parallelize_simple
#define UMAPPP_CUSTOM_PARALLEL baysor::umappp_parallel_range

#include "knncolle/knncolle.hpp"
#include "umappp/umappp.hpp"

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

namespace baysor {

// ============================================================================
// umap_embed — from raw data matrix
// ============================================================================

Eigen::MatrixXd umap_embed(
    const Eigen::MatrixXd& data,
    int ndim_out,
    int n_neighbors,
    int n_epochs,
    int seed,
    double spread,
    double min_dist
) {
    int ndim_in = static_cast<int>(data.rows());
    int nobs    = static_cast<int>(data.cols());

    if (nobs == 0 || ndim_in == 0) return Eigen::MatrixXd(ndim_out, 0);
    n_neighbors = std::min(n_neighbors, nobs - 1);

    // Build KNN via knncolle (Vantage-point tree, Euclidean).
    knncolle::SimpleMatrix<int, int, double> mat(ndim_in, nobs, data.data());
    auto index = knncolle::VptreeBuilder<knncolle::EuclideanDistance>().build_unique(mat);

    // Build neighbor list using Searcher API (knncolle v2.3+). The per-query
    // search is embarrassingly parallel: results are disjoint per index, so
    // the loop runs on the Baysor pool with one searcher per chunk.
    knncolle::NeighborList<int, double> neighbors(nobs);
    run_parallel_chunks(0, nobs, 256, Scheduling::Dynamic,
        [&](std::int64_t b, std::int64_t e, int) {
        auto searcher = index->initialize();
        std::vector<int>    out_idx;
        std::vector<double> out_dist;
        for (int i = static_cast<int>(b); i < static_cast<int>(e); ++i) {
            searcher->search(i, n_neighbors, &out_idx, &out_dist);
            neighbors[i].reserve(out_idx.size());
            for (size_t j = 0; j < out_idx.size(); ++j) {
                neighbors[i].push_back({out_idx[j], out_dist[j]});
            }
        }
    });

    // Random initialization of the embedding (RANDOM init to avoid irlba).
    std::vector<double> emb_buf(ndim_out * nobs);
    {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> ud(-10.0, 10.0);
        for (auto& v : emb_buf) v = ud(rng);
    }

    umappp::Options opt;
    opt.num_epochs = n_epochs;
    opt.seed       = static_cast<uint64_t>(seed);
    opt.spread     = spread;
    opt.min_dist   = min_dist;
    opt.num_threads = thread_pool_size();
    // Keep umappp's parallel layout optimizer off: its output is deterministic
    // and matches the serial optimizer at any thread count, but it spawns its
    // own busy-wait worker threads, which measured 1.8x slower at 4 threads
    // and 4.3x slower at 8 on a busy host. The KNN search and similarity
    // smoothing run on the Baysor pool through UMAPPP_CUSTOM_PARALLEL.
    opt.parallel_optimization = false;
    opt.initialize = umappp::InitializeMethod::NONE; // use our pre-filled buffer

    auto status = umappp::initialize(std::move(neighbors), ndim_out, emb_buf.data(), opt);
    status.run();

    // Copy result into Eigen matrix (ndim_out x nobs, column-major).
    return Eigen::Map<Eigen::MatrixXd>(emb_buf.data(), ndim_out, nobs);
}

// ============================================================================
// umap_embed_precomputed — from symmetric distance matrix
// ============================================================================

Eigen::MatrixXd umap_embed_precomputed(
    const Eigen::MatrixXd& dist_mat,
    int ndim_out,
    int n_neighbors,
    int n_epochs,
    int seed
) {
    int n = static_cast<int>(dist_mat.rows());
    if (n == 0) return Eigen::MatrixXd(ndim_out, 0);
    n_neighbors = std::min(n_neighbors, n - 1);

    // Extract k-nearest neighbors per row from the distance matrix.
    knncolle::NeighborList<int, double> neighbors(n);
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) {
        std::iota(order.begin(), order.end(), 0);
        // Partial sort to find the k smallest distances (excluding self).
        std::nth_element(order.begin(), order.begin() + n_neighbors, order.end(),
            [&](int a, int b) {
                return dist_mat(i, a) < dist_mat(i, b);
            });
        neighbors[i].reserve(n_neighbors);
        for (int j = 0; j < n_neighbors; ++j) {
            int nb = order[j];
            if (nb == i) {
                // Include the next closest instead — keep exactly n_neighbors.
                nb = order[n_neighbors]; // one past the partition boundary
            }
            neighbors[i].push_back({nb, dist_mat(i, nb)});
        }
        // Sort by distance ascending (expected by umappp).
        std::sort(neighbors[i].begin(), neighbors[i].end(),
            [](const auto& a, const auto& b) { return a.second < b.second; });
    }

    // Random initialization.
    std::vector<double> emb_buf(ndim_out * n);
    {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> ud(-10.0, 10.0);
        for (auto& v : emb_buf) v = ud(rng);
    }

    umappp::Options opt;
    opt.num_epochs = n_epochs;
    opt.seed       = static_cast<uint64_t>(seed);
    opt.min_dist   = 0.1;
    opt.spread     = 1.0;
    opt.num_threads = thread_pool_size();
    opt.parallel_optimization = false;
    opt.initialize = umappp::InitializeMethod::NONE;

    auto status = umappp::initialize(std::move(neighbors), ndim_out, emb_buf.data(), opt);
    status.run();

    return Eigen::Map<Eigen::MatrixXd>(emb_buf.data(), ndim_out, n);
}

} // namespace baysor
