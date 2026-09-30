#include "baysor/processing/data_processing/neighborhood_composition.h"
#include "baysor/utils/thread_pool.h"

#include <third_party/nanoflann.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>

namespace baysor {

namespace {

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

class KnnBlockSearcher {
public:
    explicit KnnBlockSearcher(const Eigen::MatrixXd& tree_points)
        : n_dims_(static_cast<int>(tree_points.rows()))
        , n_tree_(static_cast<int>(tree_points.cols()))
        , adaptor_(tree_points)
        , tree_(n_dims_, adaptor_, nanoflann::KDTreeSingleIndexAdaptorParams(/* max_leaf = */ 10)) {}

    template<class Callback>
    void for_each_block(
        const Eigen::MatrixXd& query_points,
        int k,
        bool sorted,
        int block_size,
        Callback&& callback
    ) const {
        const int n_query = static_cast<int>(query_points.cols());
        if (n_tree_ == 0 || n_query == 0 || k <= 0) return;

        k = std::min(k, n_tree_);
        for (int block_start = 0; block_start < n_query; block_start += block_size) {
            const int block_n = std::min(block_size, n_query - block_start);
            std::vector<int> indices(static_cast<size_t>(block_n) * static_cast<size_t>(k));
            std::vector<double> distances(static_cast<size_t>(block_n) * static_cast<size_t>(k));

            parallel_for(0, block_n, 256, [&](int local_i) {
                const int global_i = block_start + local_i;
                int* idx_ptr = indices.data() + static_cast<size_t>(local_i) * static_cast<size_t>(k);
                double* dist_ptr = distances.data() + static_cast<size_t>(local_i) * static_cast<size_t>(k);

                nanoflann::KNNResultSet<double, int> result_set(k);
                result_set.init(idx_ptr, dist_ptr);
                tree_.findNeighbors(
                    result_set,
                    query_points.col(global_i).data(),
                    nanoflann::SearchParameters(/* eps = */ 0.0f, /* sorted = */ sorted)
                );

                if (sorted) {
                    std::vector<int> order(k);
                    std::iota(order.begin(), order.end(), 0);
                    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
                        if (dist_ptr[a] != dist_ptr[b]) return dist_ptr[a] < dist_ptr[b];
                        return idx_ptr[a] < idx_ptr[b];
                    });

                    std::vector<int> sorted_indices(k);
                    std::vector<double> sorted_distances(k);
                    for (int j = 0; j < k; ++j) {
                        sorted_indices[j] = idx_ptr[order[j]];
                        sorted_distances[j] = dist_ptr[order[j]];
                    }
                    std::copy(sorted_indices.begin(), sorted_indices.end(), idx_ptr);
                    std::copy(sorted_distances.begin(), sorted_distances.end(), dist_ptr);
                }

                for (int j = 0; j < k; ++j) {
                    dist_ptr[j] = std::sqrt(dist_ptr[j]);
                }
            });

            callback(block_start, block_n, k, indices, distances);
        }
    }

    template<class Callback>
    void for_each_index_block(
        const Eigen::MatrixXd& query_points,
        const std::vector<int>& query_ids,
        int k,
        bool sorted,
        int block_size,
        Callback&& callback
    ) const {
        const int n_query = static_cast<int>(query_ids.size());
        if (n_tree_ == 0 || n_query == 0 || k <= 0) return;

        k = std::min(k, n_tree_);
        for (int block_start = 0; block_start < n_query; block_start += block_size) {
            const int block_n = std::min(block_size, n_query - block_start);
            std::vector<int> indices(static_cast<size_t>(block_n) * static_cast<size_t>(k));
            std::vector<double> distances(static_cast<size_t>(block_n) * static_cast<size_t>(k));

            parallel_for(0, block_n, 256, [&](int local_i) {
                const int query_id = query_ids[block_start + local_i];
                int* idx_ptr = indices.data() + static_cast<size_t>(local_i) * static_cast<size_t>(k);
                double* dist_ptr = distances.data() + static_cast<size_t>(local_i) * static_cast<size_t>(k);

                nanoflann::KNNResultSet<double, int> result_set(k);
                result_set.init(idx_ptr, dist_ptr);
                tree_.findNeighbors(
                    result_set,
                    query_points.col(query_id).data(),
                    nanoflann::SearchParameters(/* eps = */ 0.0f, /* sorted = */ sorted)
                );

                if (sorted) {
                    std::vector<int> order(k);
                    std::iota(order.begin(), order.end(), 0);
                    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
                        if (dist_ptr[a] != dist_ptr[b]) return dist_ptr[a] < dist_ptr[b];
                        return idx_ptr[a] < idx_ptr[b];
                    });

                    std::vector<int> sorted_indices(k);
                    std::vector<double> sorted_distances(k);
                    for (int j = 0; j < k; ++j) {
                        sorted_indices[j] = idx_ptr[order[j]];
                        sorted_distances[j] = dist_ptr[order[j]];
                    }
                    std::copy(sorted_indices.begin(), sorted_indices.end(), idx_ptr);
                    std::copy(sorted_distances.begin(), sorted_distances.end(), dist_ptr);
                }

                for (int j = 0; j < k; ++j) {
                    dist_ptr[j] = std::sqrt(dist_ptr[j]);
                }
            });

            callback(block_start, block_n, k, indices, distances);
        }
    }

private:
    int n_dims_;
    int n_tree_;
    EigenColMajorAdaptor adaptor_;
    KDTree tree_;
};

struct NeighborhoodScratch {
    std::vector<int> marks;
    std::vector<int> bucket_marks;
    std::vector<int> touched_buckets;
    std::vector<float> counts;
    int current_mark = 1;
    static constexpr int bucket_width = 64;

    void ensure(int n_genes) {
        if (static_cast<int>(marks.size()) != n_genes) {
            marks.assign(n_genes, 0);
            counts.assign(n_genes, 0.0f);
            bucket_marks.assign((n_genes + bucket_width - 1) / bucket_width, 0);
            touched_buckets.clear();
            current_mark = 1;
        }
    }

    int next_mark() {
        if (current_mark == std::numeric_limits<int>::max()) {
            std::fill(marks.begin(), marks.end(), 0); // GCOVR_EXCL_LINE: impractical to test: needs >= 2^31 queries on one thread
            current_mark = 1; // GCOVR_EXCL_LINE: impractical to test: needs >= 2^31 queries on one thread
        }
        return current_mark++;
    }
};

std::vector<int> make_query_ids(int n, const std::vector<int>* query_ids) {
    if (query_ids) return *query_ids;
    std::vector<int> ids(n);
    std::iota(ids.begin(), ids.end(), 0);
    return ids;
}

std::vector<double> exact_closest_nonzero_distances(
    const KnnBlockSearcher& searcher,
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& query_ids,
    int block_size
) {
    const int n_query = static_cast<int>(query_ids.size());
    const int n_total = static_cast<int>(pos_data.cols());
    std::vector<double> closest_nonzero(n_query, 1e-15);
    std::vector<int> unresolved(n_query);
    std::iota(unresolved.begin(), unresolved.end(), 0);

    int current_k = std::min(2, n_total);
    while (!unresolved.empty()) {
        std::vector<int> next_unresolved;
        next_unresolved.reserve(unresolved.size() / 8 + 1);

        std::vector<int> unresolved_query_ids(unresolved.size());
        for (size_t i = 0; i < unresolved.size(); ++i) unresolved_query_ids[i] = query_ids[unresolved[i]];

        searcher.for_each_index_block(pos_data, unresolved_query_ids, current_k, true, block_size,
            [&](int block_start, int block_n, int block_k,
                const std::vector<int>& indices,
                const std::vector<double>& distances) {
                // Per-index results are written to disjoint slots; the
                // "unresolved" set is collected through a flag array and
                // rebuilt in index order below, so the outcome does not
                // depend on scheduling or on the thread count.
                std::vector<unsigned char> is_unresolved(block_n, 0);

                parallel_for(0, block_n, 256, [&](int local_i) {
                    const int global_i = unresolved[block_start + local_i];
                    const double* dist_ptr = distances.data() + static_cast<size_t>(local_i) * static_cast<size_t>(block_k);

                    double closest = 1e-15;
                    for (int j = 0; j < block_k; ++j) {
                        if (dist_ptr[j] > 1e-15) {
                            closest = dist_ptr[j];
                            break;
                        }
                    }

                    if (closest > 1e-15 || current_k == n_total) {
                        closest_nonzero[global_i] = closest;
                    } else {
                        is_unresolved[local_i] = 1;
                    }
                });

                std::vector<int> local_unresolved;
                local_unresolved.reserve(static_cast<size_t>(block_n));
                for (int local_i = 0; local_i < block_n; ++local_i) {
                    if (is_unresolved[local_i]) {
                        local_unresolved.push_back(unresolved[block_start + local_i]);
                    }
                }

                next_unresolved.insert(next_unresolved.end(),
                                       local_unresolved.begin(), local_unresolved.end());
            });

        if (next_unresolved.empty() || current_k == n_total) break;
        unresolved.swap(next_unresolved);
        current_k = std::min(n_total, current_k * 2);
    }

    return closest_nonzero;
}

template<class Callback>
void for_each_neighborhood_block_with_searcher(
    const KnnBlockSearcher& searcher,
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const std::vector<int>& query_ids,
    int k,
    int n_genes,
    const std::vector<double>* confidences,
    bool normalize_by_dist,
    bool normalize,
    double distance_floor,
    Callback&& callback
) {
    constexpr int block_size = 32768;
    double med_closest_dist = distance_floor;
    if (normalize_by_dist && med_closest_dist <= 0.0) {
        std::vector<double> closest_nonzero =
            exact_closest_nonzero_distances(searcher, pos_data, query_ids, block_size);
        std::vector<double> positives;
        positives.reserve(closest_nonzero.size());
        for (double d : closest_nonzero) {
            if (d > 1e-15) positives.push_back(d);
        }
        med_closest_dist = 1e-15;
        if (!positives.empty()) {
            auto mid = positives.begin() + static_cast<std::ptrdiff_t>(positives.size() / 2);
            std::nth_element(positives.begin(), mid, positives.end());
            med_closest_dist = *mid;
        }
    }

    searcher.for_each_index_block(pos_data, query_ids, k, true, block_size,
        [&](int block_start, int block_n, int block_k,
            const std::vector<int>& indices,
            const std::vector<double>& distances) {
            using StorageIndex = Eigen::SparseMatrix<float>::StorageIndex;
            std::vector<std::vector<StorageIndex>> block_rows(block_n);
            std::vector<std::vector<float>> block_vals(block_n);

            parallel_for(0, block_n, 256, [&](int local_i) {
                thread_local NeighborhoodScratch scratch;
                scratch.ensure(n_genes);
                const int mark = scratch.next_mark();
                scratch.touched_buckets.clear();
                int nnz = 0;

                    const int* idx_ptr = indices.data() + static_cast<size_t>(local_i) * static_cast<size_t>(block_k);
                    const double* dist_ptr = distances.data() + static_cast<size_t>(local_i) * static_cast<size_t>(block_k);

                    for (int j = 0; j < block_k; ++j) {
                        const int idx = idx_ptr[j];
                        const int gene = genes[idx] - 1;
                        if (gene < 0 || gene >= n_genes) continue;

                        if (scratch.marks[gene] != mark) {
                            scratch.marks[gene] = mark;
                            scratch.counts[gene] = 0.0f;
                            ++nnz;
                            const int bucket = gene / NeighborhoodScratch::bucket_width;
                            if (scratch.bucket_marks[bucket] != mark) {
                                scratch.bucket_marks[bucket] = mark;
                                scratch.touched_buckets.push_back(bucket);
                            }
                        }

                        float weight = 1.0f;
                        if (normalize_by_dist) {
                            weight = static_cast<float>(1.0 / std::max(dist_ptr[j], med_closest_dist));
                        } else if (confidences) {
                            weight = static_cast<float>((*confidences)[idx]);
                        }
                        scratch.counts[gene] += weight;
                    }

                    float total = 0.0f;
                    if (normalize) {
                        for (int bucket : scratch.touched_buckets) {
                            const int g0 = bucket * NeighborhoodScratch::bucket_width;
                            const int g1 = std::min(n_genes, g0 + NeighborhoodScratch::bucket_width);
                            for (int gene = g0; gene < g1; ++gene) {
                                if (scratch.marks[gene] == mark) total += scratch.counts[gene];
                            }
                        }
                    }

                    std::sort(scratch.touched_buckets.begin(), scratch.touched_buckets.end());
                    auto& rows = block_rows[local_i];
                    auto& vals = block_vals[local_i];
                    rows.reserve(static_cast<size_t>(nnz));
                    vals.reserve(static_cast<size_t>(nnz));

                    for (int bucket : scratch.touched_buckets) {
                        const int g0 = bucket * NeighborhoodScratch::bucket_width;
                        const int g1 = std::min(n_genes, g0 + NeighborhoodScratch::bucket_width);
                        for (int gene = g0; gene < g1; ++gene) {
                            if (scratch.marks[gene] != mark) continue;
                            float value = scratch.counts[gene];
                            scratch.counts[gene] = 0.0f;
                            if (normalize && total > 0.0f) value /= total;
                            if (value > 1e-5f) {
                                rows.push_back(static_cast<StorageIndex>(gene));
                                vals.push_back(value);
                            }
                        }
                    }
            });

            callback(block_start, block_rows, block_vals);
        });
}

template<class Callback>
void for_each_neighborhood_block(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const std::vector<int>& query_ids,
    int k,
    int n_genes,
    const std::vector<double>* confidences,
    bool normalize_by_dist,
    bool normalize,
    double distance_floor,
    Callback&& callback
) {
    KnnBlockSearcher searcher(pos_data);
    for_each_neighborhood_block_with_searcher(
        searcher, pos_data, genes, query_ids, k, n_genes, confidences,
        normalize_by_dist, normalize, distance_floor,
        std::forward<Callback>(callback)
    );
}

Eigen::MatrixXf dense_times_sparse(
    const Eigen::MatrixXf& left,
    const Eigen::SparseMatrix<float>& right
) {
    const int n_components = static_cast<int>(left.rows());
    const int n_cols = static_cast<int>(right.cols());
    Eigen::MatrixXf out = Eigen::MatrixXf::Zero(n_components, n_cols);

    parallel_for(0, n_cols, 256, [&](int col) {
        auto out_col = out.col(col);
        for (Eigen::SparseMatrix<float>::InnerIterator it(right, col); it; ++it) {
            out_col.noalias() += it.value() * left.col(it.row());
        }
    });
    return out;
}

Eigen::MatrixXf dense_times_sparse_transpose(
    const Eigen::MatrixXf& left,
    const Eigen::SparseMatrix<float>& right
) {
    const int n_components = static_cast<int>(left.rows());
    const int n_rows = static_cast<int>(right.rows());
    const int n_cols = static_cast<int>(right.cols());

    // Deterministic chunked reduction: per-bucket accumulators merged in index
    // order. With 1 thread this is a single sequential accumulation, exactly
    // matching the former serial path.
    auto combine = [](const Eigen::MatrixXf& a, const Eigen::MatrixXf& b) {
        Eigen::MatrixXf r = a;
        r += b;
        return r;
    };
    Eigen::MatrixXf out = parallel_reduce<Eigen::MatrixXf>(
        0, n_cols, /*bucket_size=*/256, Eigen::MatrixXf::Zero(n_components, n_rows),
        [&](std::int64_t b, std::int64_t e, Eigen::MatrixXf& local) {
            for (std::int64_t col = b; col < e; ++col) {
                auto left_col = left.col(col);
                for (Eigen::SparseMatrix<float>::InnerIterator it(right, static_cast<int>(col)); it; ++it) {
                    local.col(it.row()).noalias() += it.value() * left_col;
                }
            }
        },
        combine
    );
    return out;
}

void compute_diag_and_total_var(
    const Eigen::SparseMatrix<float>& count_matrix,
    Eigen::VectorXf& diag_vals,
    Eigen::VectorXf& total_var
) {
    const int n_genes = static_cast<int>(count_matrix.rows());
    const int n_mols = static_cast<int>(count_matrix.cols());

    std::vector<float> col_sums(n_mols, 0.0f);
    parallel_for(0, n_mols, 256, [&](int col) {
        float sum = 0.0f;
        for (Eigen::SparseMatrix<float>::InnerIterator it(count_matrix, col); it; ++it) {
            sum += it.value();
        }
        col_sums[col] = sum;
    });

    // Deterministic chunked reduction over (diag, total) accumulators; see
    // dense_times_sparse_transpose for the rationale.
    using Acc = std::pair<Eigen::VectorXf, Eigen::VectorXf>;
    auto combine = [](const Acc& a, const Acc& b) {
        Acc r = a;
        r.first += b.first;
        r.second += b.second;
        return r;
    };
    Acc totals = parallel_reduce<Acc>(
        0, n_mols, /*bucket_size=*/256,
        Acc{Eigen::VectorXf::Zero(n_genes), Eigen::VectorXf::Zero(n_genes)},
        [&](std::int64_t b, std::int64_t e, Acc& acc) {
            auto& diag_local = acc.first;
            auto& total_local = acc.second;
            for (std::int64_t col = b; col < e; ++col) {
                const float col_sum = col_sums[col];
                for (Eigen::SparseMatrix<float>::InnerIterator it(count_matrix, static_cast<int>(col)); it; ++it) {
                    const int row = it.row();
                    const float value = it.value();
                    diag_local(row) += value * value;
                    total_local(row) += value * col_sum;
                }
            }
        },
        combine
    );

    diag_vals = std::move(totals.first);
    total_var = std::move(totals.second);
}

} // namespace

// ============================================================================
// neighborhood_count_matrix
// ============================================================================

Eigen::SparseMatrix<float> neighborhood_count_matrix(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int k,
    int n_genes,
    const std::vector<double>* confidences,
    bool normalize_by_dist,
    bool normalize
) {
    int n = static_cast<int>(pos_data.cols());
    if (n == 0) return Eigen::SparseMatrix<float>(n_genes, 0);
    std::vector<int> query_ids = make_query_ids(n, nullptr);
    return neighborhood_count_matrix_subset(
        pos_data, genes, query_ids, k, n_genes, confidences, normalize_by_dist, normalize
    );
}

Eigen::SparseMatrix<float> neighborhood_count_matrix_subset(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const std::vector<int>& query_ids,
    int k,
    int n_genes,
    const std::vector<double>* confidences,
    bool normalize_by_dist,
    bool normalize,
    double distance_floor
) {
    const int n_query = static_cast<int>(query_ids.size());
    const int n_total = static_cast<int>(pos_data.cols());
    if (n_query == 0) return Eigen::SparseMatrix<float>(n_genes, 0);

    if (k < 3) k = 3;
    if (!normalize) normalize_by_dist = false;
    if (confidences && normalize_by_dist) confidences = nullptr;

    k = std::min(k, n_total);
    if (n_genes <= 0) n_genes = *std::max_element(genes.begin(), genes.end());

    using StorageIndex = Eigen::SparseMatrix<float>::StorageIndex;
    std::vector<StorageIndex> outer(static_cast<size_t>(n_query) + 1, 0);
    std::vector<StorageIndex> inner;
    std::vector<float> values;

    for_each_neighborhood_block(
        pos_data, genes, query_ids, k, n_genes, confidences,
        normalize_by_dist, normalize, distance_floor,
        [&](int block_start, const auto& block_rows, const auto& block_vals) {
            size_t block_nnz = 0;
            for (const auto& rows : block_rows) block_nnz += rows.size();
            inner.reserve(inner.size() + block_nnz);
            values.reserve(values.size() + block_nnz);
            for (size_t local_i = 0; local_i < block_rows.size(); ++local_i) {
                const size_t query_idx = static_cast<size_t>(block_start) + local_i;
                outer[query_idx + 1] =
                    outer[query_idx] + static_cast<StorageIndex>(block_rows[local_i].size());
                inner.insert(inner.end(), block_rows[local_i].begin(), block_rows[local_i].end());
                values.insert(values.end(), block_vals[local_i].begin(), block_vals[local_i].end());
            }
        });

    Eigen::Map<const Eigen::SparseMatrix<float>> mapped(
        n_genes, n_query,
        static_cast<Eigen::Index>(values.size()),
        outer.data(), inner.data(), values.data());
    Eigen::SparseMatrix<float> result = mapped;
    return result;
}

double neighborhood_distance_floor(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>* query_ids
) {
    const int n = static_cast<int>(pos_data.cols());
    if (n == 0) return 1e-15;
    constexpr int block_size = 32768;
    KnnBlockSearcher searcher(pos_data);
    std::vector<int> ids = make_query_ids(n, query_ids);
    std::vector<double> closest_nonzero =
        exact_closest_nonzero_distances(searcher, pos_data, ids, block_size);
    std::vector<double> positives;
    positives.reserve(closest_nonzero.size());
    for (double d : closest_nonzero) if (d > 1e-15) positives.push_back(d);
    if (positives.empty()) return 1e-15;
    auto mid = positives.begin() + static_cast<std::ptrdiff_t>(positives.size() / 2);
    std::nth_element(positives.begin(), mid, positives.end());
    return *mid;
}

// ============================================================================
// estimate_gene_vectors (randomized indexing)
// ============================================================================

Eigen::MatrixXf estimate_gene_vectors(
    const Eigen::SparseMatrix<float>& count_matrix,
    const std::vector<int>& gene_ids,
    int n_components,
    const std::string& method,
    bool per_molecule
) {
    int n_genes = static_cast<int>(count_matrix.rows());
    int n_mols = static_cast<int>(count_matrix.cols());
    (void)gene_ids;
    (void)method;

    if (n_genes == 0 || n_mols == 0) return Eigen::MatrixXf();

    // Random projection matrix, stored as components x genes so dense_times_sparse
    // can work column-wise against the sparse count matrix.
    std::mt19937 rng(42);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    Eigen::MatrixXf random_vectors_t(n_components, n_genes);
    for (int i = 0; i < n_genes; ++i)
        for (int j = 0; j < n_components; ++j)
            random_vectors_t(j, i) = normal(rng);

    // Avoid explicitly forming coexpr = count_matrix * count_matrix.transpose().
    // We only need:
    //   1. diag(coexpr)
    //   2. row_sums(coexpr)
    //   3. coexpr * random_vectors
    // All of those can be derived by streaming over sparse columns.
    Eigen::VectorXf diag_vals;
    Eigen::VectorXf total_var;
    compute_diag_and_total_var(count_matrix, diag_vals, total_var);

    constexpr float var_clip = 0.05f;
    Eigen::VectorXf clipped_diag = diag_vals;
    if (var_clip > 0 && n_genes > 1) {
        std::vector<float> df_sorted(n_genes);
        for (int i = 0; i < n_genes; ++i)
            df_sorted[i] = (total_var(i) > 0) ? diag_vals(i) / total_var(i) : 0.0f;
        std::sort(df_sorted.begin(), df_sorted.end());
        float q = df_sorted[static_cast<int>((1.0f - var_clip) * (n_genes - 1))];

        for (int i = 0; i < n_genes; ++i)
            clipped_diag(i) = std::min(q * total_var(i), diag_vals(i));
    }

    Eigen::VectorXf delta = diag_vals - clipped_diag;
    Eigen::VectorXf row_sums = total_var - delta;

    // Unclipped coexpr * random_vectors = count_matrix * (count_matrix^T * random_vectors).
    Eigen::MatrixXf proj_mol = dense_times_sparse(random_vectors_t, count_matrix); // k x n_mols
    Eigen::MatrixXf gene_emb_t = dense_times_sparse_transpose(proj_mol, count_matrix); // k x n_genes

    // Diagonal clipping only changes the diagonal of coexpr, so subtract the removed
    // diagonal mass directly from each gene vector, then normalize by the clipped row sums.
    for (int i = 0; i < n_genes; ++i) {
        gene_emb_t.col(i).noalias() -= delta(i) * random_vectors_t.col(i);
        if (row_sums(i) > 0) gene_emb_t.col(i) /= row_sums(i);
    }

    if (per_molecule) {
        return dense_times_sparse(gene_emb_t, count_matrix);
    }

    return gene_emb_t;  // n_components x n_genes
}

Eigen::MatrixXf project_gene_vectors(
    const Eigen::MatrixXf& gene_emb_t,
    const Eigen::SparseMatrix<float>& count_matrix
) {
    return dense_times_sparse(gene_emb_t, count_matrix);
}

Eigen::MatrixXf project_neighborhood_vectors(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int k,
    const Eigen::MatrixXf& gene_emb_t,
    int n_genes,
    const std::vector<int>* query_ids,
    const std::vector<double>* confidences,
    bool normalize_by_dist,
    bool normalize,
    double distance_floor,
    bool log_transform
) {
    const int n_total = static_cast<int>(pos_data.cols());
    std::vector<int> ids = make_query_ids(n_total, query_ids);
    if (ids.empty()) return Eigen::MatrixXf(gene_emb_t.rows(), 0);
    Eigen::MatrixXf out = Eigen::MatrixXf::Zero(gene_emb_t.rows(), static_cast<int>(ids.size()));
    stream_projected_neighborhood_vectors(
        pos_data, genes, k, gene_emb_t, n_genes, &ids, confidences,
        normalize_by_dist, normalize, distance_floor, log_transform, 32768,
        [&](int block_start, const std::vector<int>&, const Eigen::MatrixXf& block_vecs) {
            out.middleCols(block_start, block_vecs.cols()) = block_vecs;
        }
    );
    return out;
}

void stream_projected_neighborhood_vectors(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int k,
    const Eigen::MatrixXf& gene_emb_t,
    int n_genes,
    const std::vector<int>* query_ids,
    const std::vector<double>* confidences,
    bool normalize_by_dist,
    bool normalize,
    double distance_floor,
    bool log_transform,
    int block_size,
    const std::function<void(int, const std::vector<int>&, const Eigen::MatrixXf&)>& callback
) {
    const int n_total = static_cast<int>(pos_data.cols());
    std::vector<int> ids = make_query_ids(n_total, query_ids);
    if (ids.empty() || !callback) return;
    if (n_genes <= 0) n_genes = *std::max_element(genes.begin(), genes.end());
    if (confidences && normalize_by_dist) confidences = nullptr;
    k = std::min(std::max(k, 3), n_total);

    KnnBlockSearcher searcher(pos_data);
    double med_closest_dist = distance_floor;
    if (normalize_by_dist && med_closest_dist <= 0.0) {
        std::vector<double> closest_nonzero =
            exact_closest_nonzero_distances(searcher, pos_data, ids, block_size);
        std::vector<double> positives;
        positives.reserve(closest_nonzero.size());
        for (double d : closest_nonzero) if (d > 1e-15) positives.push_back(d);
        med_closest_dist = 1e-15;
        if (!positives.empty()) {
            auto mid = positives.begin() + static_cast<std::ptrdiff_t>(positives.size() / 2);
            std::nth_element(positives.begin(), mid, positives.end());
            med_closest_dist = *mid;
        }
    }

    searcher.for_each_index_block(pos_data, ids, k, true, block_size,
        [&](int block_start, int block_n, int block_k,
            const std::vector<int>& indices,
            const std::vector<double>& distances) {
            Eigen::MatrixXf block_vecs = Eigen::MatrixXf::Zero(gene_emb_t.rows(), block_n);

            parallel_for(0, block_n, 256, [&](int local_i) {
                thread_local NeighborhoodScratch scratch;
                scratch.ensure(n_genes);
                const int mark = scratch.next_mark();
                scratch.touched_buckets.clear();

                const int* idx_ptr = indices.data() + static_cast<size_t>(local_i) * static_cast<size_t>(block_k);
                    const double* dist_ptr = distances.data() + static_cast<size_t>(local_i) * static_cast<size_t>(block_k);

                    for (int j = 0; j < block_k; ++j) {
                        const int idx = idx_ptr[j];
                        const int gene = genes[idx] - 1;
                        if (gene < 0 || gene >= n_genes) continue;
                        if (scratch.marks[gene] != mark) {
                            scratch.marks[gene] = mark;
                            scratch.counts[gene] = 0.0f;
                            const int bucket = gene / NeighborhoodScratch::bucket_width;
                            if (scratch.bucket_marks[bucket] != mark) {
                                scratch.bucket_marks[bucket] = mark;
                                scratch.touched_buckets.push_back(bucket);
                            }
                        }

                        float weight = 1.0f;
                        if (normalize_by_dist) {
                            weight = static_cast<float>(1.0 / std::max(dist_ptr[j], med_closest_dist));
                        } else if (confidences) {
                            weight = static_cast<float>((*confidences)[idx]);
                        }
                        scratch.counts[gene] += weight;
                    }

                    float total = 0.0f;
                    if (normalize) {
                        for (int bucket : scratch.touched_buckets) {
                            const int g0 = bucket * NeighborhoodScratch::bucket_width;
                            const int g1 = std::min(n_genes, g0 + NeighborhoodScratch::bucket_width);
                            for (int gene = g0; gene < g1; ++gene) {
                                if (scratch.marks[gene] == mark) total += scratch.counts[gene];
                            }
                        }
                    }

                    std::sort(scratch.touched_buckets.begin(), scratch.touched_buckets.end());
                    Eigen::VectorXf col = Eigen::VectorXf::Zero(gene_emb_t.rows());
                    for (int bucket : scratch.touched_buckets) {
                        const int g0 = bucket * NeighborhoodScratch::bucket_width;
                        const int g1 = std::min(n_genes, g0 + NeighborhoodScratch::bucket_width);
                        for (int gene = g0; gene < g1; ++gene) {
                            if (scratch.marks[gene] != mark) continue;
                            float value = scratch.counts[gene];
                            scratch.counts[gene] = 0.0f;
                            if (normalize && total > 0.0f) value /= total;
                            if (log_transform) value = static_cast<float>(std::log(value * 10000.0f + 1e-5f));
                            if (value > 1e-5f) col.noalias() += value * gene_emb_t.col(gene);
                        }
                    }
                    block_vecs.col(local_i) = col;
            });

            std::vector<int> block_query_ids(
                ids.begin() + block_start, // GCOVR_EXCL_LINE: dead GCC block; expression counted on line 809
                ids.begin() + block_start + block_n
            );
            callback(block_start, block_query_ids, block_vecs);
        });
}

} // namespace baysor
