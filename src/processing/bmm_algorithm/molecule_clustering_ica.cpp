// ICA initialisation of the MRF molecule clustering (cluster_molecules_ica):
// gene co-occurrence matrix -> whitening -> FastICA -> initial expression
// profiles. Kept apart from the MRF EM (molecule_clustering.cpp) so that the
// irlba templates do not change how the EM loop is compiled.

#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/utils/thread_pool.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <irlba/compute.hpp>
#include <irlba/wrappers.hpp>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace baysor {

// ============================================================================
// ICA whitening: top-k eigenpairs of the row covariance
//   C = Xc * Xc' / (n - 1),  Xc = X - rowmean(X)
// ============================================================================

namespace detail {

// Dense path, matching Julia's MultivariateStats.ICA: full covariance and a
// complete symmetric eigen-decomposition, of which the top k pairs are kept.
// O(n_features^2 * n_samples + n_features^3) time, five dense
// n_features x n_features matrices.
IcaWhitening ica_whitening_dense(const Eigen::MatrixXd& X, int k) {
    const int n_features = static_cast<int>(X.rows());
    const int n_samples  = static_cast<int>(X.cols());
    Eigen::VectorXd mean_vec = X.rowwise().mean();
    Eigen::MatrixXd Xc = X.colwise() - mean_vec;
    Eigen::MatrixXd cov =
        (Xc * Xc.transpose()) / static_cast<double>(std::max(1, n_samples - 1));
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
    IcaWhitening out;
    out.eigenvalues.resize(k);
    out.eigenvectors.resize(n_features, k);
    for (int i = 0; i < k; ++i) {
        int src_idx = n_features - 1 - i;
        out.eigenvalues(i) = eig.eigenvalues()(src_idx);
        out.eigenvectors.col(i) = eig.eigenvectors().col(src_idx);
    }
    return out;
}

} // namespace detail

namespace {

// irlba matrix interface (irlba/MockMatrix.hpp) for A = X' with X sparse
// (CSC, plus its transpose so that both products are column dot products).
// irlba centres the columns of A, i.e. the rows of X, so the right singular
// vectors of the centred A are the eigenvectors of C. Every output entry is
// summed in a fixed order by one task, so the products do not depend on the
// thread count.
class SparseTransposeOperator {
public:
    using SpMat = Eigen::SparseMatrix<double>;

    SparseTransposeOperator(const SpMat& x, const SpMat& x_t) : x_(x), x_t_(x_t) {}

    Eigen::Index rows() const { return x_.cols(); }
    Eigen::Index cols() const { return x_.rows(); }

    struct Workspace {};
    struct AdjointWorkspace {};
    Workspace workspace() const { return {}; }
    AdjointWorkspace adjoint_workspace() const { return {}; }

    // out = X' * rhs
    template <class Right>
    void multiply(const Right& rhs, Workspace&, Eigen::VectorXd& out) const {
        column_dots(x_, rhs, out);
    }

    // out = X * rhs
    template <class Right>
    void adjoint_multiply(const Right& rhs, AdjointWorkspace&, Eigen::VectorXd& out) const {
        column_dots(x_t_, rhs, out);
    }

    template <class EigenMatrix>
    EigenMatrix realize() const { return EigenMatrix(x_t_); }

private:
    // out(j) = sum over the stored entries (i, j) of m of m(i, j) * rhs(i)
    template <class Right>
    static void column_dots(const SpMat& m, const Right& rhs, Eigen::VectorXd& out) {
        out.resize(m.cols());
        parallel_for(0, m.cols(), 64, [&](std::int64_t j) {
            double sum = 0.0;
            for (SpMat::InnerIterator it(m, static_cast<Eigen::Index>(j)); it; ++it)
                sum += it.value() * rhs.coeff(it.row());
            out(j) = sum;
        });
    }

    const SpMat& x_;
    const SpMat& x_t_;
};

// Row means of X from the CSC of its transpose.
Eigen::VectorXd sparse_row_means(const Eigen::SparseMatrix<double>& x_t) {
    Eigen::VectorXd means(x_t.cols());
    const double n = static_cast<double>(std::max<Eigen::Index>(1, x_t.rows()));
    for (Eigen::Index i = 0; i < x_t.cols(); ++i) {
        double sum = 0.0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(x_t, i); it; ++it) sum += it.value();
        means(i) = sum / n;
    }
    return means;
}

// Top-k eigenpairs of C from a truncated SVD (irlba) of the implicitly
// centred A = X' (rows = samples): eigenvalues d^2 / (n - 1), eigenvectors
// the right singular vectors. Each vector's sign is fixed so that its
// largest-magnitude entry is positive (the first one on ties), which makes
// the FastICA start independent of the solver's sign choice.
template <class Operator>
detail::IcaWhitening truncated_whitening(const Operator& op, const Eigen::VectorXd& row_means, int k) {
    const int n_samples = static_cast<int>(op.rows());
    irlba::Centered<Operator, Eigen::VectorXd> centered(op, row_means);
    irlba::Options opt;
    opt.convergence_tolerance = 1e-10;
    opt.max_iterations = 10000;
    opt.seed = 42;
    Eigen::MatrixXd U, V;
    Eigen::VectorXd D;
    const auto stats = irlba::compute(centered, k, U, V, D, opt);
    if (!stats.first)
        throw std::runtime_error("truncated eigen-decomposition for the ICA whitening did not converge");
    spdlog::debug("ICA whitening: truncated eigen-decomposition converged after {} restarts.", stats.second);

    const double denom = static_cast<double>(std::max(1, n_samples - 1));
    detail::IcaWhitening out;
    out.eigenvalues.resize(k);
    out.eigenvectors.resize(op.cols(), k);
    for (int i = 0; i < k; ++i) {
        out.eigenvalues(i) = D(i) * D(i) / denom;
        out.eigenvectors.col(i) = V.col(i);
        Eigen::Index imax = 0;
        out.eigenvectors.col(i).cwiseAbs().maxCoeff(&imax);
        if (out.eigenvectors(imax, i) < 0) out.eigenvectors.col(i) *= -1.0;
    }
    return out;
}

// W0 = P * diag(1 / sqrt(lambda)) with tiny eigenvalues clamped to avoid
// division by ~0.
Eigen::MatrixXd whitening_matrix(detail::IcaWhitening wh) {
    for (Eigen::Index i = 0; i < wh.eigenvalues.size(); ++i)
        if (wh.eigenvalues(i) < 1e-10) wh.eigenvalues(i) = 1e-10;
    return wh.eigenvectors * wh.eigenvalues.cwiseSqrt().cwiseInverse().asDiagonal(); // m × k
}

// Symmetric FastICA (tanh nonlinearity) on whitened data Z (k × n_samples).
// Returns the rotation W (k × k).
Eigen::MatrixXd fast_ica_rotation(
    const Eigen::MatrixXd& Z,
    int n_components,
    int max_iter,
    double tol,
    unsigned int seed
) {
    const int n_samples = static_cast<int>(Z.cols());
    std::mt19937 rng(seed);
    std::normal_distribution<double> ndist(0.0, 1.0);
    Eigen::MatrixXd W(n_components, n_components);
    for (int c = 0; c < n_components; ++c) {
        for (int r = 0; r < n_components; ++r) {
            W(r, c) = ndist(rng);
        }
        double norm = W.col(c).norm();
        if (norm > 0.0) W.col(c) /= norm;
    }

    Eigen::MatrixXd W_prev(n_components, n_components);
    Eigen::MatrixXd U(n_samples, n_components);
    Eigen::MatrixXd Y(n_components, n_components);
    Eigen::VectorXd E1(n_components);

    for (int iter = 0; iter < max_iter; ++iter) {
        W_prev = W;

        // U <- W' * X, stored as (n_samples × n_components)
        U.noalias() = Z.transpose() * W;

        // Tanh nonlinearity with a=1 and mean derivative per component.
        for (int c = 0; c < n_components; ++c) {
            double deriv_sum = 0.0;
            for (int i = 0; i < n_samples; ++i) {
                double t = std::tanh(U(i, c));
                U(i, c) = t;
                deriv_sum += 1.0 - t * t;
            }
            E1(c) = deriv_sum / static_cast<double>(n_samples);
        }

        // Y <- E{x g(w'x)}
        Y.noalias() = (Z * U) / static_cast<double>(n_samples);

        // W <- Y - E{g'(w'x)} * W
        for (int c = 0; c < n_components; ++c) {
            W.col(c) = Y.col(c) - E1(c) * W.col(c);
        }

        // Symmetric decorrelation: W <- W * (W'W)^(-1/2)
        Eigen::MatrixXd gram = W.transpose() * W;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> gram_eig(gram);
        Eigen::VectorXd gram_vals = gram_eig.eigenvalues();
        for (int i = 0; i < gram_vals.size(); ++i) {
            if (gram_vals(i) < 1e-10) gram_vals(i) = 1e-10;
        }
        Eigen::MatrixXd invsqrt =
            gram_eig.eigenvectors() *
            gram_vals.cwiseSqrt().cwiseInverse().asDiagonal() *
            gram_eig.eigenvectors().transpose();
        W = W * invsqrt;

        Eigen::MatrixXd prod = W * W_prev.transpose();
        double chg = 0.0;
        for (int i = 0; i < n_components; ++i) {
            chg = std::max(chg, std::abs(std::abs(prod(i, i)) - 1.0));
        }
        if (chg < tol) break;
    }
    return W;
}

// Julia's fit(ICA, X, k) rejects k > min(m, n) with an error that
// cluster_molecules_on_mrf's wrapper catches to fall back to hash/random
// init. Clamping instead would return fewer columns than the caller expects
// and make it index W out of bounds.
void check_ica_components(Eigen::Index n_features, Eigen::Index n_samples, int n_components) {
    if (n_components > std::min(n_features, n_samples))
        throw std::invalid_argument("k must not exceed min(m, n).");
}

} // namespace

namespace detail {

IcaWhitening ica_whitening_truncated(const Eigen::SparseMatrix<double>& X, int k) {
    check_ica_components(X.rows(), X.cols(), k);
    const Eigen::SparseMatrix<double> X_t = X.transpose();
    return truncated_whitening(SparseTransposeOperator(X, X_t), sparse_row_means(X_t), k);
}

Eigen::SparseMatrix<double> sparse_gene_spatial_cor(
    const std::vector<int>& genes,
    const std::vector<double>& confidence,
    const AdjList& adj_list,
    double confidence_threshold
) {
    using SpMat = Eigen::SparseMatrix<double>;
    const int n = static_cast<int>(genes.size());
    const int n_genes = *std::max_element(genes.begin(), genes.end());  // 1-based max

    // Molecules above the confidence threshold, grouped by gene in molecule order.
    std::vector<int> gene_offsets(static_cast<size_t>(n_genes) + 1, 0);
    for (int gi = 0; gi < n; ++gi) {
        if (confidence[gi] < confidence_threshold) continue;
        const int g = genes[gi] - 1;
        if (g >= 0 && g < n_genes) ++gene_offsets[static_cast<size_t>(g) + 1];
    }
    std::partial_sum(gene_offsets.begin(), gene_offsets.end(), gene_offsets.begin());
    std::vector<int> mols_by_gene(static_cast<size_t>(gene_offsets.back()));
    {
        std::vector<int> next(gene_offsets.begin(), gene_offsets.end() - 1);
        for (int gi = 0; gi < n; ++gi) {
            if (confidence[gi] < confidence_threshold) continue;
            const int g = genes[gi] - 1;
            if (g >= 0 && g < n_genes) mols_by_gene[static_cast<size_t>(next[g]++)] = gi;
        }
    }

    // Row g2 (the column g2 of the transpose): summed neighbour edge weights
    // per neighbour gene, accumulated in the same order as the dense builder,
    // into a per-worker dense scratch row; only the touched genes are stored.
    std::vector<std::vector<int>> row_genes(static_cast<size_t>(n_genes));
    std::vector<std::vector<double>> row_vals(static_cast<size_t>(n_genes));
    const int n_workers = thread_pool_size();
    std::vector<std::vector<double>> scratch(static_cast<size_t>(n_workers));
    std::vector<std::vector<int>> marks(static_cast<size_t>(n_workers));
    parallel_for(0, n_genes, 1, [&](int g2, int w) {
        auto& row = scratch[static_cast<size_t>(w)];
        auto& mark = marks[static_cast<size_t>(w)];
        if (row.empty()) {
            row.assign(static_cast<size_t>(n_genes), 0.0);
            mark.assign(static_cast<size_t>(n_genes), -1);
        }
        auto& touched = row_genes[static_cast<size_t>(g2)];
        for (int j = gene_offsets[g2]; j < gene_offsets[g2 + 1]; ++j) {
            const int gi = mols_by_gene[static_cast<size_t>(j)];
            const int nc = adj_list.neighbor_count(gi);
            const int32_t* nb_ids = adj_list.neighbor_ids(gi);
            const double* nb_wts = adj_list.neighbor_weights(gi);
            for (int ai = 0; ai < nc; ++ai) {
                const int nb = nb_ids[ai];
                if (confidence[nb] < confidence_threshold) continue;
                const int g1 = genes[nb] - 1;
                if (g1 < 0 || g1 >= n_genes) continue;
                if (mark[static_cast<size_t>(g1)] != g2) {
                    mark[static_cast<size_t>(g1)] = g2;
                    row[static_cast<size_t>(g1)] = 0.0;
                    touched.push_back(g1);
                }
                row[static_cast<size_t>(g1)] += nb_wts[ai];
            }
        }
        std::sort(touched.begin(), touched.end());
        auto& vals = row_vals[static_cast<size_t>(g2)];
        vals.reserve(touched.size());
        for (int g1 : touched) vals.push_back(row[static_cast<size_t>(g1)]);
    });

    // CSC of the transpose (column g2 = row g2), then the matrix itself.
    SpMat x_t(n_genes, n_genes);
    {
        size_t nnz = 0;
        for (const auto& r : row_genes) nnz += r.size();
        x_t.resizeNonZeros(static_cast<Eigen::Index>(nnz));
        auto* outer = x_t.outerIndexPtr();
        Eigen::Index pos = 0;
        for (int g2 = 0; g2 < n_genes; ++g2) {
            const auto& r = row_genes[static_cast<size_t>(g2)];
            const auto& v = row_vals[static_cast<size_t>(g2)];
            std::copy(r.begin(), r.end(), x_t.innerIndexPtr() + pos);
            std::copy(v.begin(), v.end(), x_t.valuePtr() + pos);
            pos += static_cast<Eigen::Index>(r.size());
            outer[g2 + 1] = static_cast<SpMat::StorageIndex>(pos);
            std::vector<int>().swap(row_genes[static_cast<size_t>(g2)]);
            std::vector<double>().swap(row_vals[static_cast<size_t>(g2)]);
        }
    }
    SpMat x = x_t.transpose();

    // Normalise by sqrt(sum_weight[r] * sum_weight[c]), sum_weight = row sum +
    // column sum, floored at 0.1 as in the dense builder (zeros stay zero).
    std::vector<double> sum_weight(static_cast<size_t>(n_genes), 0.0);
    parallel_for(0, n_genes, 64, [&](int g) {
        double row_sum = 0.0, col_sum = 0.0;
        for (SpMat::InnerIterator it(x_t, g); it; ++it) row_sum += it.value();
        for (SpMat::InnerIterator it(x, g); it; ++it) col_sum += it.value();
        sum_weight[static_cast<size_t>(g)] = row_sum + col_sum;
    });
    parallel_for(0, n_genes, 64, [&](int c) {
        for (SpMat::InnerIterator it(x, c); it; ++it) {
            const double denom = std::sqrt(sum_weight[static_cast<size_t>(it.row())] * sum_weight[static_cast<size_t>(c)]);
            it.valueRef() /= std::max(denom, 0.1);
        }
    });
    return x;
}

} // namespace detail

// ============================================================================
// fast_ica — symmetric deflation FastICA (tanh nonlinearity)
//
// Input:  X  — data matrix (n_features × n_samples)
//         n_components — number of independent components
// Returns: unmixing matrix W (n_features × n_components)
//          whose columns correspond to the components used by Julia's
//          MultivariateStats.ICA fit via ica_fit.W.
//
// Mirrors Julia's MultivariateStats.ICA usage in cluster_molecules_on_mrf,
// including its argument check `k <= min(m, n) || error(...)`: throws
// std::invalid_argument instead of silently clamping, so the caller's
// try/catch falls back exactly like Julia's wrapper does.
//
// ============================================================================
static Eigen::MatrixXd fast_ica(
    const Eigen::MatrixXd& X,
    int n_components,
    int max_iter = 1000,
    double tol   = 1e-5,
    unsigned int seed = 42
) {
    check_ica_components(X.rows(), X.cols(), n_components);

    // 1. Center rows to match Julia's preprocess_mean/centralize path.
    // 2. Whiten with the top-k covariance eigenvectors, matching:
    //    C = Z * Z' / (n - 1); W0 = P * Diagonal(1 ./ sqrt.(v)); Z = W0' * Z
    Eigen::MatrixXd W0 = whitening_matrix(detail::ica_whitening_dense(X, n_components));
    Eigen::VectorXd mean_vec = X.rowwise().mean();
    Eigen::MatrixXd Xc = X.colwise() - mean_vec;
    Eigen::MatrixXd Z = W0.transpose() * Xc;                // k × n

    // 3. FastICA rotation; 4. map the whitened solution back to the original
    //    feature space.
    return W0 * fast_ica_rotation(Z, n_components, max_iter, tol, seed);
}

// fast_ica for a sparse X with the truncated whitening. It computes the same
// subspace as the dense eigen-decomposition; the result differs in
// floating-point rounding and in the (arbitrary) eigenvector signs, which
// change the FastICA start point. Memory and time are O(nnz * k) per Lanczos
// step instead of O(n^3) with five dense n x n matrices.
static Eigen::MatrixXd fast_ica_sparse(
    const Eigen::SparseMatrix<double>& X,
    int n_components,
    int max_iter = 1000,
    double tol   = 1e-5,
    unsigned int seed = 42
) {
    check_ica_components(X.rows(), X.cols(), n_components);
    const Eigen::SparseMatrix<double> X_t = X.transpose();
    const Eigen::VectorXd mean_vec = sparse_row_means(X_t);
    Eigen::MatrixXd W0 = whitening_matrix(
        truncated_whitening(SparseTransposeOperator(X, X_t), mean_vec, n_components));

    // Z = W0' * (X - mean 1') without forming the dense centred X.
    const Eigen::VectorXd shift = W0.transpose() * mean_vec;
    Eigen::MatrixXd Z(n_components, X.cols());
    parallel_for(0, X.cols(), 64, [&](std::int64_t j) {
        Eigen::VectorXd col = -shift;
        for (Eigen::SparseMatrix<double>::InnerIterator it(X, static_cast<Eigen::Index>(j)); it; ++it)
            col.noalias() += it.value() * W0.row(it.row()).transpose();
        Z.col(j) = col;
    });
    return W0 * fast_ica_rotation(Z, n_components, max_iter, tol, seed);
}

// ============================================================================
// cluster_molecules_ica
//
// Port of Julia's DataFrame-based wrapper for cluster_molecules_on_mrf:
//   1. Compute pairwise gene spatial co-occurrence (already in color_utils)
//   2. Run FastICA on the correlation matrix to get n_clusters gene profiles
//   3. Call core EM with ICA init; fall back to hash init if ICA throws
// ============================================================================
ClusteringResult cluster_molecules_ica(
    const std::vector<int>& genes,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    int n_clusters,
    double tol,
    double mrf_weight,
    int max_iters,
    bool verbose,
    int dense_whitening_max_genes
) {
    if (n_clusters <= 1) return {};

    // Infer n_genes
    int n_genes = 0;
    for (int g : genes) if (g > n_genes) n_genes = g;

    // 1. Gene spatial correlation matrix (n_genes × n_genes), 0-based genes.
    // 2. FastICA on the correlation matrix → unmixing matrix (n_genes × n_clusters).
    //    Large panels use a sparse matrix (gene pairs that are never
    //    neighbours stay zero; ~0.1 % of the pairs are non-zero on a
    //    whole-transcriptome panel) and the truncated whitening: the dense
    //    one is cubic in the number of genes (hours at ~10,000 genes).
    const bool truncated = n_genes > dense_whitening_max_genes;
    Eigen::MatrixXd cor_mat;
    Eigen::SparseMatrix<double> sparse_cor_mat;
    if (truncated) {
        sparse_cor_mat = detail::sparse_gene_spatial_cor(genes, confidence, adj_list);
        if (verbose)
            spdlog::info("ICA whitening: truncated eigen-decomposition ({} genes > {}, {} non-zero co-occurrences).",
                         n_genes, dense_whitening_max_genes, sparse_cor_mat.nonZeros());
    } else {
        cor_mat = pairwise_gene_spatial_cor(genes, confidence, adj_list);
    }

    std::unique_ptr<Eigen::MatrixXd> exprs_init_ptr;
    try {
        Eigen::MatrixXd W = truncated ? fast_ica_sparse(sparse_cor_mat, n_clusters)
                                      : fast_ica(cor_mat, n_clusters);
        // Convert mixing matrix to expression profiles:
        //   ct_exprs_init[k][g] = abs(W[g][k]) / sum_g'(abs(W[g'][k]))
        // Matches Julia: (abs.(ica_fit.W) ./ sum(abs.(ica_fit.W), dims=1))'
        Eigen::MatrixXd exprs(n_clusters, n_genes);
        for (int k = 0; k < n_clusters; ++k) {
            double col_sum = W.col(k).cwiseAbs().sum();
            if (col_sum < 1e-10) col_sum = 1.0;
            for (int g = 0; g < n_genes; ++g)
                exprs(k, g) = std::abs(W(g, k)) / col_sum;
        }
        exprs_init_ptr = std::make_unique<Eigen::MatrixXd>(std::move(exprs));
        if (verbose) spdlog::info("ICA initialization succeeded ({} components).", n_clusters);
    } catch (const std::exception& e) {
        // Discard any partially built init so *every* failure path falls back
        // to hash initialization, consistent with the log message (and with
        // Julia, which leaves ct_exprs_init = nothing on any exception).
        exprs_init_ptr.reset();
        spdlog::warn("ICA did not converge ({}), falling back to hash initialization.", e.what());
    } catch (...) {
        exprs_init_ptr.reset(); // GCOVR_EXCL_LINE: this handler is entered only by a non-std exception, and no Baysor or dependency code throws one
        spdlog::warn("ICA failed, falling back to hash initialization."); // GCOVR_EXCL_LINE: reachable only if a non-std exception escapes the ICA try block; same justification as main.cpp's catch-all (no Baysor or dependency exception lacks std::exception)
    }

    // 3. Core EM with ICA init (or nullptr → hash fallback)
    return cluster_molecules_on_mrf(
        genes, adj_list, confidence,
        n_clusters, tol, mrf_weight, max_iters, verbose,
        exprs_init_ptr.get()
    );
}

} // namespace baysor
