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
#include "baysor/processing/data_processing/umap_optimize.h"

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

namespace baysor {

namespace {

// The fuzzy simplicial set of a kNN list, as umappp::initialize() builds it:
// smoothed similarities, then the symmetrised union.
void knn_to_fuzzy_graph(knncolle::NeighborList<int, double>& neighbors) {
    const umappp::Options opt;
    umappp::internal::neighbor_similarities<int, double>(neighbors, opt.local_connectivity, opt.bandwidth);
    umappp::internal::combine_neighbor_sets<int, double>(neighbors, opt.mix_ratio);
}

} // namespace

// ============================================================================
// umap_fuzzy_graph / umap_embed_graph / umap_embed — from raw data matrix
// ============================================================================

UmapGraph umap_fuzzy_graph(const Eigen::MatrixXd& data, int n_neighbors) {
    int ndim_in = static_cast<int>(data.rows());
    int nobs    = static_cast<int>(data.cols());

    if (nobs == 0 || ndim_in == 0) return {};
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

    knn_to_fuzzy_graph(neighbors);
    return neighbors;
}

// umappp::initialize(graph, ...) with InitializeMethod::NONE on a uniform
// random start, followed by Status::run(), with the layout optimised by
// Baysor's serial optimiser (umap_optimize.h) instead of umappp's. The set-up
// steps are umappp's own, in umappp's order. The optimiser differs from
// umappp's serial one only in computing pow(d2, b) with fast_pow (relative
// error < 1e-13), which is much cheaper than glibc's pow and changes the NCV
// colours by mean dE ~3. umappp's parallel optimiser is not used: it is
// deterministic and matches the serial one at any thread count, but it spawns
// its own busy-wait worker threads, which measured 1.8x slower at 4 threads
// and 4.3x slower at 8 on a busy host.
Eigen::MatrixXd umap_embed_graph(
    const UmapGraph& graph,
    int ndim_out,
    int n_epochs,
    int seed,
    double spread,
    double min_dist
) {
    const int nobs = static_cast<int>(graph.size());
    if (nobs == 0) return Eigen::MatrixXd(ndim_out, 0);

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

    const auto ab = umappp::internal::find_ab(opt.spread, opt.min_dist);
    const int num_epochs = umappp::internal::choose_num_epochs(opt.num_epochs, graph.size());
    auto epochs = umappp::internal::similarities_to_epochs<int, double>(
        graph, num_epochs, opt.negative_sample_rate);
    std::mt19937_64 engine(opt.seed);
    umap_detail::optimize_layout_dispatch</*FastPow_=*/true>(
        static_cast<std::size_t>(ndim_out), emb_buf.data(), epochs,
        ab.first, ab.second, opt.repulsion_strength, opt.learning_rate,
        engine, epochs.total_epochs);

    // Copy result into Eigen matrix (ndim_out x nobs, column-major).
    return Eigen::Map<Eigen::MatrixXd>(emb_buf.data(), ndim_out, nobs);
}

Eigen::MatrixXd umap_embed(
    const Eigen::MatrixXd& data,
    int ndim_out,
    int n_neighbors,
    int n_epochs,
    int seed,
    double spread,
    double min_dist
) {
    return umap_embed_graph(umap_fuzzy_graph(data, n_neighbors), ndim_out, n_epochs, seed, spread, min_dist);
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

    knn_to_fuzzy_graph(neighbors);
    return umap_embed_graph(neighbors, ndim_out, n_epochs, seed, /*spread=*/1.0, /*min_dist=*/0.1);
}

} // namespace baysor
