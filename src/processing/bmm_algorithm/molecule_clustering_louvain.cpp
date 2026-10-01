#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/processing/data_processing/neighborhood_composition.h"
#include "baysor/utils/thread_pool.h"

#include <spdlog/spdlog.h>
#include <third_party/nanoflann.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

namespace baysor {

namespace {

struct LouvainGraph {
    AdjList adj;
    std::vector<double> degree;
    double total_weight_twice = 0.0;
};

struct GraphStats {
    int n_nodes = 0;
    int n_isolated = 0;
    int n_components = 0;
    int largest_component = 0;
};

struct PartitionAttempt {
    double resolution = 1.0;
    std::vector<int> membership;
    std::vector<double> move_fracs;
    int n_clusters = 0;
};

static LouvainGraph make_louvain_graph(const AdjList& graph) {
    LouvainGraph out;
    out.adj = graph;
    out.degree.resize(graph.n_molecules(), 0.0);
    for (int i = 0; i < graph.n_molecules(); ++i) {
        double deg = 0.0;
        const int nc = graph.neighbor_count(i);
        const double* wt = graph.neighbor_weights(i);
        for (int j = 0; j < nc; ++j) deg += wt[j];
        out.degree[i] = deg;
        out.total_weight_twice += deg;
    }
    return out;
}

static GraphStats compute_graph_stats(const AdjList& graph) {
    GraphStats stats;
    stats.n_nodes = graph.n_molecules();
    if (stats.n_nodes == 0) return stats;

    std::vector<unsigned char> seen(stats.n_nodes, 0);
    std::vector<int> stack;
    stack.reserve(1024);

    for (int i = 0; i < stats.n_nodes; ++i) {
        const int degree = graph.neighbor_count(i);
        if (degree == 0) {
            ++stats.n_isolated;
            if (!seen[i]) {
                seen[i] = 1;
                ++stats.n_components;
                stats.largest_component = std::max(stats.largest_component, 1);
            }
            continue;
        }
        if (seen[i]) continue;

        int comp_size = 0;
        ++stats.n_components;
        stack.clear();
        stack.push_back(i);
        seen[i] = 1;
        while (!stack.empty()) {
            const int cur = stack.back();
            stack.pop_back();
            ++comp_size;
            const int nc = graph.neighbor_count(cur);
            const int32_t* nb_ids = graph.neighbor_ids(cur);
            for (int j = 0; j < nc; ++j) {
                const int nb = nb_ids[j];
                if (!seen[nb]) {
                    seen[nb] = 1;
                    stack.push_back(nb);
                }
            }
        }
        stats.largest_component = std::max(stats.largest_component, comp_size);
    }

    return stats;
}

static std::vector<int> reindex_membership_zero_based(const std::vector<int>& membership, int* n_clusters_out = nullptr) {
    std::unordered_map<int, int> remap;
    remap.reserve(membership.size());
    std::vector<int> out(membership.size(), 0);
    int next_id = 0;
    for (size_t i = 0; i < membership.size(); ++i) {
        auto it = remap.find(membership[i]);
        if (it == remap.end()) {
            it = remap.emplace(membership[i], next_id++).first;
        }
        out[i] = it->second;
    }
    if (n_clusters_out) *n_clusters_out = next_id;
    return out;
}

static double community_edge_gain(
    int node,
    int target_comm,
    double resolution,
    const std::vector<double>& comm_tot,
    const std::vector<double>& neigh_w,
    const LouvainGraph& graph
) {
    return neigh_w[target_comm] -
           resolution * graph.degree[node] * comm_tot[target_comm] /
               std::max(graph.total_weight_twice, 1e-12);
}

static bool louvain_one_level(
    const LouvainGraph& graph,
    double resolution,
    std::vector<int>& membership,
    std::vector<double>* move_fracs = nullptr,
    int max_passes = 100
) {
    const int n = graph.adj.n_molecules();
    if (n == 0) return false;

    int n_comms = 0;
    membership = reindex_membership_zero_based(membership, &n_comms);
    std::vector<double> comm_tot(n_comms, 0.0);
    for (int i = 0; i < n; ++i) comm_tot[membership[i]] += graph.degree[i];

    bool moved_any = false;
    std::vector<double> neigh_w(std::max(n_comms, 1), 0.0);
    std::vector<int> touched;
    touched.reserve(64);

    for (int pass = 0; pass < max_passes; ++pass) {
        int moved = 0;

        for (int node = 0; node < n; ++node) {
            const int cur_comm = membership[node];
            const int nc = graph.adj.neighbor_count(node);
            const int32_t* nb_ids = graph.adj.neighbor_ids(node);
            const double* nb_wts = graph.adj.neighbor_weights(node);

            touched.clear();
            for (int ai = 0; ai < nc; ++ai) {
                const int nb = nb_ids[ai];
                if (nb == node) continue; // self-loops are constant w.r.t. moves
                const int comm = membership[nb];
                if (neigh_w[comm] == 0.0) touched.push_back(comm);
                neigh_w[comm] += nb_wts[ai];
            }

            comm_tot[cur_comm] -= graph.degree[node];

            int best_comm = cur_comm;
            double best_gain = 0.0;

            for (int comm : touched) {
                double gain = community_edge_gain(node, comm, resolution, comm_tot, neigh_w, graph);
                if (gain > best_gain + 1e-12 ||
                    (std::abs(gain - best_gain) <= 1e-12 && comm < best_comm)) {
                    best_gain = gain;
                    best_comm = comm;
                }
            }

            membership[node] = best_comm;
            comm_tot[best_comm] += graph.degree[node];
            if (best_comm != cur_comm) {
                moved_any = true;
                ++moved;
            }

            for (int comm : touched) neigh_w[comm] = 0.0;
        }

        if (move_fracs) move_fracs->push_back(static_cast<double>(moved) / std::max(n, 1));
        if (moved == 0) break;
    }

    return moved_any;
}

static std::vector<int> split_disconnected_communities_zero_based(
    const AdjList& graph,
    const std::vector<int>& membership
) {
    const int n = graph.n_molecules();
    if (n == 0) return {};

    std::vector<int> refined(n, -1);
    std::vector<unsigned char> seen(n, 0);
    std::vector<int> stack;
    stack.reserve(256);
    int next_comm = 0;

    for (int i = 0; i < n; ++i) {
        if (seen[i]) continue;
        const int orig_comm = membership[i];
        seen[i] = 1;
        stack.clear();
        stack.push_back(i);
        while (!stack.empty()) {
            const int cur = stack.back();
            stack.pop_back();
            refined[cur] = next_comm;
            const int nc = graph.neighbor_count(cur);
            const int32_t* nb_ids = graph.neighbor_ids(cur);
            for (int j = 0; j < nc; ++j) {
                const int nb = nb_ids[j];
                if (seen[nb] || membership[nb] != orig_comm) continue;
                seen[nb] = 1;
                stack.push_back(nb);
            }
        }
        ++next_comm;
    }

    return refined;
}

static AdjList aggregate_graph(
    const LouvainGraph& graph,
    const std::vector<int>& membership,
    int n_comms
) {
    std::unordered_map<std::uint64_t, double> edge_weight;
    edge_weight.reserve(static_cast<size_t>(graph.adj.nnz()));

    auto key = [](int a, int b) -> std::uint64_t {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
               static_cast<std::uint32_t>(b);
    };

    for (int i = 0; i < graph.adj.n_molecules(); ++i) {
        const int ci = membership[i];
        const int nc = graph.adj.neighbor_count(i);
        const int32_t* nb_ids = graph.adj.neighbor_ids(i);
        const double* nb_wts = graph.adj.neighbor_weights(i);

        for (int ai = 0; ai < nc; ++ai) {
            const int j = nb_ids[ai];
            const double w = nb_wts[ai];
            if (j < i) continue;

            const int cj = membership[j];
            if (i == j) {
                edge_weight[key(ci, ci)] += w;
            } else if (ci == cj) {
                edge_weight[key(ci, ci)] += 2.0 * w;
            } else {
                const int a = std::min(ci, cj);
                const int b = std::max(ci, cj);
                edge_weight[key(a, b)] += w;
            }
        }
    }

    AdjList agg;
    agg.indptr.assign(n_comms + 1, 0);

    int nnz = 0;
    for (const auto& kv : edge_weight) {
        int a = static_cast<int>(kv.first >> 32);
        int b = static_cast<int>(kv.first & 0xffffffffu);
        if (a == b) {
            agg.indptr[a + 1] += 1;
            ++nnz;
        } else {
            agg.indptr[a + 1] += 1;
            agg.indptr[b + 1] += 1;
            nnz += 2;
        }
    }

    for (int i = 1; i <= n_comms; ++i) agg.indptr[i] += agg.indptr[i - 1];
    agg.indices.resize(nnz);
    agg.weights.resize(nnz);

    std::vector<int32_t> write_pos = agg.indptr;
    for (const auto& kv : edge_weight) {
        int a = static_cast<int>(kv.first >> 32);
        int b = static_cast<int>(kv.first & 0xffffffffu);
        double w = kv.second;
        int pos = write_pos[a]++;
        agg.indices[pos] = b;
        agg.weights[pos] = w;
        if (a != b) {
            pos = write_pos[b]++;
            agg.indices[pos] = a;
            agg.weights[pos] = w;
        }
    }

    return agg;
}

static std::vector<int> run_louvain_zero_based(
    const AdjList& input_graph,
    double resolution,
    int max_passes,
    std::vector<double>* move_fracs = nullptr
) {
    LouvainGraph graph = make_louvain_graph(input_graph);
    const int n0 = graph.adj.n_molecules();
    std::vector<int> membership(n0);
    std::iota(membership.begin(), membership.end(), 0);
    std::vector<int> finest_to_current = membership;

    for (int level = 0; level < max_passes; ++level) {
        std::vector<double> level_move_fracs;
        bool moved = louvain_one_level(graph, resolution, membership,
                                       move_fracs ? &level_move_fracs : nullptr,
                                       max_passes);
        int n_comms = 0;
        membership = reindex_membership_zero_based(membership, &n_comms);

        for (int& id : finest_to_current) id = membership[id];
        if (move_fracs) {
            move_fracs->insert(move_fracs->end(), level_move_fracs.begin(), level_move_fracs.end());
        }

        if (!moved || n_comms == graph.adj.n_molecules()) break;
        graph = make_louvain_graph(aggregate_graph(graph, membership, n_comms));
        membership.resize(graph.adj.n_molecules());
        std::iota(membership.begin(), membership.end(), 0);
    }

    return reindex_membership_zero_based(finest_to_current);
}

static std::vector<int> run_leiden_zero_based(
    const AdjList& input_graph,
    double resolution,
    int max_passes,
    std::vector<double>* move_fracs = nullptr
) {
    LouvainGraph graph = make_louvain_graph(input_graph);
    const int n0 = graph.adj.n_molecules();
    std::vector<int> membership(n0);
    std::iota(membership.begin(), membership.end(), 0);
    std::vector<int> finest_to_current = membership;

    for (int level = 0; level < max_passes; ++level) {
        std::vector<double> level_move_fracs;
        bool moved = louvain_one_level(
            graph, resolution, membership,
            move_fracs ? &level_move_fracs : nullptr,
            max_passes
        );
        membership = reindex_membership_zero_based(membership);
        std::vector<int> refined = split_disconnected_communities_zero_based(graph.adj, membership);

        bool refined_changed = (refined != membership);
        int n_comms = 0;
        refined = reindex_membership_zero_based(refined, &n_comms);

        for (int& id : finest_to_current) id = refined[id];
        if (move_fracs) {
            move_fracs->insert(move_fracs->end(), level_move_fracs.begin(), level_move_fracs.end());
        }

        if ((!moved && !refined_changed) || n_comms == graph.adj.n_molecules()) break;
        graph = make_louvain_graph(aggregate_graph(graph, refined, n_comms));
        membership.resize(graph.adj.n_molecules());
        std::iota(membership.begin(), membership.end(), 0);
    }

    return reindex_membership_zero_based(finest_to_current);
}

struct EigenColMajorAdaptor {
    const Eigen::MatrixXd& mat;

    explicit EigenColMajorAdaptor(const Eigen::MatrixXd& m) : mat(m) {}

    inline size_t kdtree_get_point_count() const { return static_cast<size_t>(mat.cols()); }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return mat(static_cast<Eigen::Index>(dim), static_cast<Eigen::Index>(idx));
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, EigenColMajorAdaptor>,
    EigenColMajorAdaptor,
    -1,
    int
>;

static Eigen::MatrixXd normalize_columns_l2(const Eigen::MatrixXf& mol_vecs) {
    Eigen::MatrixXd normalized = mol_vecs.cast<double>();
    parallel_for_static(0, normalized.cols(), [&](int i) {
        const double norm = normalized.col(i).norm();
        if (norm > 1e-12) {
            normalized.col(i) /= norm;
        } else {
            normalized.col(i).setZero();
        }
    });
    return normalized;
}

} // namespace

AdjList build_knn_similarity_graph(
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int k
) {
    const int n = static_cast<int>(mol_vecs.cols());
    AdjList out;
    out.indptr.assign(n + 1, 0);
    if (n <= 1 || k <= 0) {
        return out;
    }

    k = std::min(k, n - 1);
    const int query_k = k + 1;
    Eigen::MatrixXd normalized = normalize_columns_l2(mol_vecs);

    EigenColMajorAdaptor adaptor(normalized);
    KDTree tree(
        static_cast<int>(normalized.rows()), adaptor,
        nanoflann::KDTreeSingleIndexAdaptorParams(/*max_leaf=*/10)
    );

    // Edges collected per chunk and merged in chunk order: the canonical
    // (src, neighbor) order at any thread count.
    constexpr std::int64_t chunk = 256;
    struct Edges {
        std::vector<int> src, dst;
        std::vector<double> wt;
    };
    std::vector<Edges> edges_by_chunk((n + chunk - 1) / chunk);

    run_parallel_chunks(0, n, chunk, Scheduling::Dynamic,
        [&](std::int64_t b, std::int64_t e, int) {
        auto& [src, dst, wt] = edges_by_chunk[b / chunk];
        std::vector<int> nn_indices(query_k);
        std::vector<double> nn_distances(query_k);

        for (int i = static_cast<int>(b); i < static_cast<int>(e); ++i) {
            nanoflann::KNNResultSet<double, int> result_set(query_k);
            result_set.init(nn_indices.data(), nn_distances.data());
            tree.findNeighbors(
                result_set,
                normalized.col(i).data(),
                nanoflann::SearchParameters(/*eps=*/0.0f, /*sorted=*/true)
            );

            for (int j = 0; j < query_k; ++j) {
                const int nb = nn_indices[j];
                if (nb == i) continue;
                const double cosine = std::max(0.0, normalized.col(i).dot(normalized.col(nb)));
                if (cosine <= 0.0) continue;
                const double conf_scale =
                    std::sqrt(std::max(0.0, confidence[i]) * std::max(0.0, confidence[nb]));
                const double weight = cosine * conf_scale;
                if (weight <= 0.0) continue;
                src.push_back(i);
                dst.push_back(nb);
                wt.push_back(weight);
            }
        }
    });

    Edges all;
    for (const auto& c : edges_by_chunk) {
        all.src.insert(all.src.end(), c.src.begin(), c.src.end());
        all.dst.insert(all.dst.end(), c.dst.begin(), c.dst.end());
        all.wt.insert(all.wt.end(), c.wt.begin(), c.wt.end());
    }
    const auto& [src, dst, wt] = all;

    if (src.empty()) return out;
    return AdjList::from_edge_list(src.data(), dst.data(), wt.data(),
                                   static_cast<int>(src.size()), n);
}

static std::vector<int> transfer_labels_from_anchor_vectors_exact(
    const Eigen::MatrixXf& anchor_vecs,
    const std::vector<int>& anchor_labels,
    const Eigen::MatrixXf& mol_vecs,
    int k
) {
    const int n_anchors = static_cast<int>(anchor_vecs.cols());
    const int n = static_cast<int>(mol_vecs.cols());
    std::vector<int> out(n, 1);
    if (n_anchors == 0 || n == 0 || static_cast<int>(anchor_labels.size()) != n_anchors) return out;

    const int max_label = *std::max_element(anchor_labels.begin(), anchor_labels.end());
    const int query_k = std::max(1, std::min(k, n_anchors));
    Eigen::MatrixXd normalized_anchors = normalize_columns_l2(anchor_vecs);
    EigenColMajorAdaptor adaptor(normalized_anchors);
    KDTree tree(
        static_cast<int>(normalized_anchors.rows()), adaptor,
        nanoflann::KDTreeSingleIndexAdaptorParams(/*max_leaf=*/10)
    );

    run_parallel_chunks(0, n, 512, Scheduling::Dynamic,
        [&](std::int64_t b, std::int64_t e, int) {
        std::vector<int> nn_indices(query_k);
        std::vector<double> nn_distances(query_k);
        std::vector<double> query(static_cast<size_t>(mol_vecs.rows()), 0.0);
        std::vector<double> label_w(static_cast<size_t>(max_label + 1), 0.0);
        std::vector<int> touched;
        touched.reserve(query_k);

        for (int i = static_cast<int>(b); i < static_cast<int>(e); ++i) {
            double norm_sq = 0.0;
            for (int d = 0; d < mol_vecs.rows(); ++d) {
                const double v = static_cast<double>(mol_vecs(d, i));
                query[static_cast<size_t>(d)] = v;
                norm_sq += v * v;
            }
            const double norm = std::sqrt(norm_sq);
            if (norm > 1e-12) {
                for (double& v : query) v /= norm;
            }

            nanoflann::KNNResultSet<double, int> result_set(query_k);
            result_set.init(nn_indices.data(), nn_distances.data());
            tree.findNeighbors(
                result_set,
                query.data(),
                nanoflann::SearchParameters(/*eps=*/0.0f, /*sorted=*/true)
            );

            int best_label = anchor_labels[nn_indices[0]];
            double best_weight = -1.0;
            touched.clear();
            for (int j = 0; j < query_k; ++j) {
                const int label = anchor_labels[nn_indices[j]];
                const double weight = 1.0 / std::max(nn_distances[j], 1e-12);
                if (label_w[label] == 0.0) touched.push_back(label);
                label_w[label] += weight;
                if (label_w[label] > best_weight + 1e-12 ||
                    (std::abs(label_w[label] - best_weight) <= 1e-12 && label < best_label)) {
                    best_weight = label_w[label];
                    best_label = label;
                }
            }
            out[i] = best_label;
            for (int label : touched) label_w[label] = 0.0;
        }
    });

    return out;
}

static PartitionAttempt run_graph_partition_once(
    const AdjList& graph,
    ClusterMethod method,
    double resolution,
    int max_passes
) {
    PartitionAttempt out;
    out.resolution = resolution;
    out.membership = method == ClusterMethod::Leiden
        ? run_leiden_zero_based(graph, resolution, max_passes, &out.move_fracs)
        : run_louvain_zero_based(graph, resolution, max_passes, &out.move_fracs);
    out.membership = reindex_membership_zero_based(out.membership, &out.n_clusters);
    return out;
}

static int compute_micro_target(int target_clusters) {
    if (target_clusters <= 1) return target_clusters;
    return std::clamp(target_clusters * 4, std::max(target_clusters + 4, 12), 200);
}

static std::vector<int> merge_microclusters_to_target_zero_based(
    const std::vector<int>& micro_membership,
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int target_clusters
) {
    int micro_clusters = 0;
    std::vector<int> membership = reindex_membership_zero_based(micro_membership, &micro_clusters);
    if (target_clusters <= 0 || micro_clusters <= target_clusters) return membership;

    const int dims = static_cast<int>(mol_vecs.rows());
    Eigen::MatrixXf centroids = Eigen::MatrixXf::Zero(dims, micro_clusters);
    std::vector<double> weights(micro_clusters, 0.0);
    for (int i = 0; i < static_cast<int>(membership.size()); ++i) {
        const int c = membership[i];
        const double w = std::max(1e-6, confidence.empty() ? 1.0 : confidence[i]);
        centroids.col(c) += static_cast<float>(w) * mol_vecs.col(i);
        weights[c] += w;
    }
    for (int c = 0; c < micro_clusters; ++c) {
        if (weights[c] > 0.0) centroids.col(c) /= static_cast<float>(weights[c]);
    }

    std::vector<int> rep(micro_clusters);
    std::iota(rep.begin(), rep.end(), 0);
    std::vector<unsigned char> active(micro_clusters, 1);
    int active_count = micro_clusters;

    auto cosine = [&](int a, int b) -> double {
        const float na = centroids.col(a).norm();
        const float nb = centroids.col(b).norm();
        if (na <= 1e-12f || nb <= 1e-12f) return -1.0;
        return static_cast<double>(centroids.col(a).dot(centroids.col(b)) / (na * nb));
    };

    while (active_count > target_clusters) {
        int best_a = -1;
        int best_b = -1;
        double best_sim = -std::numeric_limits<double>::infinity();
        for (int a = 0; a < micro_clusters; ++a) {
            if (!active[a]) continue;
            for (int b = a + 1; b < micro_clusters; ++b) {
                if (!active[b]) continue;
                const double sim = cosine(a, b);
                if (sim > best_sim) {
                    best_sim = sim;
                    best_a = a;
                    best_b = b;
                }
            }
        }
        if (best_a < 0 || best_b < 0) break;

        const double wa = weights[best_a];
        const double wb = weights[best_b];
        const double wsum = wa + wb;
        if (wsum > 0.0) {
            centroids.col(best_a) =
                (static_cast<float>(wa) * centroids.col(best_a) +
                 static_cast<float>(wb) * centroids.col(best_b)) /
                static_cast<float>(wsum);
        }
        weights[best_a] = wsum;
        active[best_b] = 0;
        --active_count;

        for (int& r : rep) {
            if (r == best_b) r = best_a;
        }
    }

    int final_clusters = 0;
    std::vector<int> coarse_rep = reindex_membership_zero_based(rep, &final_clusters);
    std::vector<int> coarse_membership(membership.size(), 0);
    for (int i = 0; i < static_cast<int>(membership.size()); ++i) {
        coarse_membership[i] = coarse_rep[membership[i]];
    }
    return coarse_membership;
}

std::vector<int> graph_partition_to_target(
    const AdjList& graph,
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    ClusterMethod method,
    int target_clusters,
    double resolution_seed,
    int max_passes,
    GraphClusteringSummary* summary
) {
    if (graph.n_molecules() == 0) return {};
    if (method != ClusterMethod::Louvain && method != ClusterMethod::Leiden) {
        return {};
    }

    if (target_clusters <= 1) {
        std::vector<int> out(graph.n_molecules(), 1);
        if (summary) {
            summary->micro_clusters = 1;
            summary->final_clusters = 1;
            summary->chosen_resolution = resolution_seed;
        }
        return out;
    }

    const int micro_target = compute_micro_target(target_clusters);
    std::vector<double> candidate_resolutions;
    candidate_resolutions.reserve(13);
    double seed = std::max(1e-4, resolution_seed);
    for (int s = -6; s <= 6; ++s) {
        candidate_resolutions.push_back(seed * std::pow(2.0, static_cast<double>(s)));
    }
    std::sort(candidate_resolutions.begin(), candidate_resolutions.end());
    candidate_resolutions.erase(
        std::unique(candidate_resolutions.begin(), candidate_resolutions.end(),
                    [](double a, double b) { return std::abs(a - b) <= 1e-12; }),
        candidate_resolutions.end()
    );

    std::vector<PartitionAttempt> attempts(candidate_resolutions.size());
    parallel_for(0, static_cast<int>(candidate_resolutions.size()), 1, [&](int i) {
        attempts[i] = run_graph_partition_once(
            graph, method, candidate_resolutions[static_cast<size_t>(i)], max_passes);
    });

    const PartitionAttempt* chosen = nullptr;
    for (const auto& attempt : attempts) {
        if (attempt.n_clusters >= micro_target) {
            if (chosen == nullptr || attempt.n_clusters < chosen->n_clusters) {
                chosen = &attempt;
            }
        }
    }
    if (chosen == nullptr) {
        for (const auto& attempt : attempts) {
            if (attempt.n_clusters >= target_clusters) {
                if (chosen == nullptr || attempt.n_clusters < chosen->n_clusters) {
                    chosen = &attempt;
                }
            }
        }
    }
    if (chosen == nullptr) {
        chosen = &attempts.front();
        for (const auto& attempt : attempts) {
            if (attempt.n_clusters > chosen->n_clusters) chosen = &attempt;
        }
    }

    std::vector<int> merged = merge_microclusters_to_target_zero_based(
        chosen->membership, mol_vecs, confidence, target_clusters
    );
    int final_clusters = 0;
    merged = reindex_membership_zero_based(merged, &final_clusters);
    for (int& m : merged) ++m;

    if (summary) {
        summary->micro_clusters = chosen->n_clusters;
        summary->final_clusters = final_clusters;
        summary->chosen_resolution = chosen->resolution;
        summary->move_fracs = chosen->move_fracs;
    }
    return merged;
}

std::vector<int> louvain_partition(
    const AdjList& graph,
    double resolution,
    int max_passes
) {
    std::vector<int> membership = run_louvain_zero_based(graph, resolution, max_passes, nullptr);
    for (int& m : membership) ++m;
    return membership;
}

std::vector<int> leiden_partition(
    const AdjList& graph,
    double resolution,
    int max_passes
) {
    std::vector<int> membership = run_leiden_zero_based(graph, resolution, max_passes, nullptr);
    for (int& m : membership) ++m;
    return membership;
}

static ClusteringResult cluster_molecules_graph_backend(
    ClusterMethod method,
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const std::vector<double>& confidence,
    double resolution,
    int graph_k,
    int spatial_k,
    int target_clusters,
    int n_dims,
    int basis_sample_size,
    bool verbose
) {
    if (pos_data.cols() == 0) return {};

    const int n_genes = *std::max_element(genes.begin(), genes.end());
    const int effective_spatial_k = spatial_k > 0 ? spatial_k : std::max(3, n_genes / 10);
    auto ncv_model = std::make_shared<NcvProjectedModel>(fit_ncv_projected_model(
        pos_data, genes, n_genes, confidence, effective_spatial_k, basis_sample_size, n_dims, true
    ));
    const std::vector<int>& cluster_anchor_ids = ncv_model->basis.basis_ids;
    std::vector<double> cluster_anchor_confidence;
    cluster_anchor_confidence.reserve(cluster_anchor_ids.size());
    for (int id : cluster_anchor_ids) {
        cluster_anchor_confidence.push_back(confidence.empty() ? 1.0 : confidence[id]);
    }
    const char* label = method == ClusterMethod::Leiden ? "Leiden" : "Louvain";
    if (verbose) {
        spdlog::info("{} clustering: using {} basis anchors (spatial_k={}, graph_k={}).", // GCOVR_EXCL_LINE: gcov exception-cleanup artifact
                     label, cluster_anchor_ids.size(), effective_spatial_k, graph_k);
    }
    const auto& anchor_vecs = ncv_model->basis.basis_vecs;
    AdjList weighted_graph = build_knn_similarity_graph(anchor_vecs, cluster_anchor_confidence, graph_k);
    const GraphStats graph_stats = compute_graph_stats(weighted_graph);

    GraphClusteringSummary summary;
    std::vector<int> anchor_labels = graph_partition_to_target(
        weighted_graph, anchor_vecs, cluster_anchor_confidence, method,
        target_clusters, resolution, 100, &summary
    );
    std::vector<int> membership = transfer_labels_from_anchor_vectors_exact(
        anchor_vecs, anchor_labels, ncv_model->mol_vecs, graph_k
    );

    if (verbose) {
        if (graph_stats.n_components > 1 || graph_stats.n_isolated > 0) {
            spdlog::warn(
                "{} anchor graph has {} connected components ({} isolated anchors; largest component={}).",
                label, graph_stats.n_components, graph_stats.n_isolated, graph_stats.largest_component
            );
        }
        spdlog::info(
            "{} clustering complete: {} anchor communities merged to {} clusters.",
            label, summary.micro_clusters, summary.final_clusters
        );
    }

    return ClusteringResult{
        Eigen::MatrixXd(),
        std::move(membership),
        Eigen::MatrixXd(),
        std::move(summary.move_fracs),
        {},
        std::move(ncv_model)
    };
}

ClusteringResult cluster_molecules_louvain(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const AdjList&,
    const std::vector<double>& confidence,
    double resolution,
    int graph_k,
    int spatial_k,
    int target_clusters,
    int n_dims,
    int basis_sample_size,
    bool verbose
) {
    return cluster_molecules_graph_backend(
        ClusterMethod::Louvain, pos_data, genes, confidence,
        resolution, graph_k, spatial_k, target_clusters, n_dims, basis_sample_size, verbose
    );
}

ClusteringResult cluster_molecules_leiden(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const AdjList&,
    const std::vector<double>& confidence,
    double resolution,
    int graph_k,
    int spatial_k,
    int target_clusters,
    int n_dims,
    int basis_sample_size,
    bool verbose
) {
    return cluster_molecules_graph_backend(
        ClusterMethod::Leiden, pos_data, genes, confidence,
        resolution, graph_k, spatial_k, target_clusters, n_dims, basis_sample_size, verbose
    );
}

ClusteringResult cluster_molecules(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    const ClusteringOptions& options,
    bool verbose
) {
    switch (options.method) {
        case ClusterMethod::None:
            return {};
        case ClusterMethod::Mrf:
            if (options.n_clusters <= 1) return {};
            return cluster_molecules_ica(
                genes, adj_list, confidence,
                options.n_clusters, options.tol, options.mrf_weight,
                options.max_iters, verbose
            );
        case ClusterMethod::Louvain:
            return cluster_molecules_louvain(
                pos_data, genes, adj_list, confidence,
                options.resolution, options.graph_k, options.spatial_k, options.n_clusters,
                options.n_dims, options.basis_sample_size, verbose
            );
        case ClusterMethod::Leiden:
            return cluster_molecules_leiden(
                pos_data, genes, adj_list, confidence,
                options.resolution, options.graph_k, options.spatial_k, options.n_clusters,
                options.n_dims, options.basis_sample_size, verbose
            );
    }
    return {};
}

} // namespace baysor
