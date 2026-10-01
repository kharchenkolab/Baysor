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

namespace {

using SpMat = Eigen::SparseMatrix<double>;

// irlba matrix interface (irlba/MockMatrix.hpp) for A = X' with X sparse
// (CSC, plus its transpose so that both products are column dot products).
// irlba centres the columns of A, i.e. the rows of X, so the right singular
// vectors of the centred A are the eigenvectors of the row covariance of X.
// Every output entry is summed in a fixed order by one task, so the products
// do not depend on the thread count.
class SparseTransposeOperator {
public:
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

// Row means of X, each summed in column order.
Eigen::VectorXd row_means(const SpMat& x) {
    Eigen::VectorXd sums = Eigen::VectorXd::Zero(x.rows());
    for (Eigen::Index j = 0; j < x.cols(); ++j)
        for (SpMat::InnerIterator it(x, j); it; ++it) sums(it.row()) += it.value();
    return sums / static_cast<double>(std::max<Eigen::Index>(1, x.cols()));
}

// Julia's fit(ICA, X, k) rejects k > min(m, n) with an error, on which the
// clustering falls back to the hash initialisation.
void check_ica_components(Eigen::Index n_features, Eigen::Index n_samples, int n_components) {
    if (n_components > std::min(n_features, n_samples))
        throw std::invalid_argument("k must not exceed min(m, n).");
}

// W0 = P * diag(1 / sqrt(lambda)), with tiny eigenvalues clamped.
Eigen::MatrixXd whitening_matrix(detail::IcaWhitening wh) {
    for (Eigen::Index i = 0; i < wh.eigenvalues.size(); ++i)
        if (wh.eigenvalues(i) < 1e-10) wh.eigenvalues(i) = 1e-10;
    return wh.eigenvectors * wh.eigenvalues.cwiseSqrt().cwiseInverse().asDiagonal(); // m × k
}

// Symmetric FastICA (tanh nonlinearity) on whitened data Z (k × n_samples),
// as in Julia's MultivariateStats.ICA. Returns the rotation W (k × k).
Eigen::MatrixXd fast_ica_rotation(const Eigen::MatrixXd& Z) {
    constexpr int max_iter = 1000;
    constexpr double tol = 1e-5;
    const int n_components = static_cast<int>(Z.rows());
    const int n_samples = static_cast<int>(Z.cols());
    std::mt19937 rng(42);
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

} // namespace

namespace detail {

// Julia's MultivariateStats.ICA whitening: full covariance
// C = Xc * Xc' / (n - 1), Xc = X - rowmean(X), and its complete symmetric
// eigen-decomposition, of which the top k pairs are kept.
IcaWhitening ica_whitening_dense(const Eigen::MatrixXd& X, int k) {
    check_ica_components(X.rows(), X.cols(), k);
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

// Truncated SVD of the implicitly centred A = X': eigenvalues d^2 / (n - 1),
// eigenvectors the right singular vectors. The sign of each vector is fixed
// (largest-magnitude entry positive, the first one on ties) so that the
// FastICA start does not depend on the solver's sign choice.
IcaWhitening ica_whitening_truncated(const SpMat& X, int k) {
    check_ica_components(X.rows(), X.cols(), k);
    const SpMat X_t = X.transpose();
    const SparseTransposeOperator op(X, X_t);
    const Eigen::VectorXd means = row_means(X);
    irlba::Centered<SparseTransposeOperator, Eigen::VectorXd> centered(op, means);
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

    const double denom = static_cast<double>(std::max<Eigen::Index>(1, X.cols() - 1));
    IcaWhitening out{D.head(k).cwiseAbs2() / denom, V.leftCols(k)};
    for (int i = 0; i < k; ++i) {
        Eigen::Index imax = 0;
        out.eigenvectors.col(i).cwiseAbs().maxCoeff(&imax);
        if (out.eigenvectors(imax, i) < 0) out.eigenvectors.col(i) *= -1.0;
    }
    return out;
}

SpMat sparse_gene_spatial_cor(
    const std::vector<int>& genes,
    const std::vector<double>& confidence,
    const AdjList& adj_list,
    double confidence_threshold
) {
    const int n = static_cast<int>(genes.size());
    const int n_genes = *std::max_element(genes.begin(), genes.end());

    // Molecules above the confidence threshold, grouped by gene in molecule order.
    auto skipped = [&](int i) { return confidence[i] < confidence_threshold || genes[i] <= 0; };
    std::vector<int> gene_offsets(n_genes + 1, 0);
    for (int gi = 0; gi < n; ++gi) {
        if (!skipped(gi)) ++gene_offsets[genes[gi]];
    }
    std::partial_sum(gene_offsets.begin(), gene_offsets.end(), gene_offsets.begin());
    std::vector<int> mols_by_gene(gene_offsets.back());
    {
        std::vector<int> next(gene_offsets.begin(), gene_offsets.end() - 1);
        for (int gi = 0; gi < n; ++gi) {
            if (!skipped(gi)) mols_by_gene[next[genes[gi] - 1]++] = gi;
        }
    }

    // Row g2 (the column g2 of the transpose): summed neighbour edge weights
    // per neighbour gene, accumulated in the same order as the dense builder,
    // into a per-worker dense scratch row; only the touched genes are stored.
    std::vector<std::vector<int>> row_genes(n_genes);
    std::vector<std::vector<double>> row_vals(n_genes);
    std::vector<std::vector<double>> scratch(thread_pool_size());
    std::vector<std::vector<int>> marks(thread_pool_size());
    parallel_for(0, n_genes, 1, [&](int g2, int w) {
        auto& row = scratch[w];
        auto& mark = marks[w];
        if (row.empty()) {
            row.assign(n_genes, 0.0);
            mark.assign(n_genes, -1);
        }
        auto& touched = row_genes[g2];
        for (int j = gene_offsets[g2]; j < gene_offsets[g2 + 1]; ++j) {
            const int gi = mols_by_gene[j];
            const int nc = adj_list.neighbor_count(gi);
            const int32_t* nb_ids = adj_list.neighbor_ids(gi);
            const double* nb_wts = adj_list.neighbor_weights(gi);
            for (int ai = 0; ai < nc; ++ai) {
                const int nb = nb_ids[ai];
                if (skipped(nb)) continue;
                const int g1 = genes[nb] - 1;
                if (mark[g1] != g2) {
                    mark[g1] = g2;
                    row[g1] = 0.0;
                    touched.push_back(g1);
                }
                row[g1] += nb_wts[ai];
            }
        }
        std::sort(touched.begin(), touched.end());
        auto& vals = row_vals[g2];
        vals.reserve(touched.size());
        for (int g1 : touched) vals.push_back(row[g1]);
    });

    // CSC of the transpose (column g2 = row g2), then the matrix itself.
    SpMat x_t(n_genes, n_genes);
    {
        size_t nnz = 0;
        for (const auto& r : row_genes) nnz += r.size();
        x_t.resizeNonZeros(static_cast<Eigen::Index>(nnz));
        Eigen::Index pos = 0;
        for (int g2 = 0; g2 < n_genes; ++g2) {
            std::copy(row_genes[g2].begin(), row_genes[g2].end(), x_t.innerIndexPtr() + pos);
            std::copy(row_vals[g2].begin(), row_vals[g2].end(), x_t.valuePtr() + pos);
            pos += static_cast<Eigen::Index>(row_genes[g2].size());
            x_t.outerIndexPtr()[g2 + 1] = static_cast<SpMat::StorageIndex>(pos);
            std::vector<int>().swap(row_genes[g2]);
            std::vector<double>().swap(row_vals[g2]);
        }
    }
    SpMat x = x_t.transpose();

    // Normalise by sqrt(sum_weight[r] * sum_weight[c]), sum_weight = row sum +
    // column sum, floored at 0.1 as in the dense builder (zeros stay zero).
    std::vector<double> sum_weight(n_genes, 0.0);
    parallel_for(0, n_genes, 64, [&](int g) {
        double row_sum = 0.0, col_sum = 0.0;
        for (SpMat::InnerIterator it(x_t, g); it; ++it) row_sum += it.value();
        for (SpMat::InnerIterator it(x, g); it; ++it) col_sum += it.value();
        sum_weight[g] = row_sum + col_sum;
    });
    parallel_for(0, n_genes, 64, [&](int c) {
        for (SpMat::InnerIterator it(x, c); it; ++it)
            it.valueRef() /= std::max(std::sqrt(sum_weight[it.row()] * sum_weight[c]), 0.1);
    });
    return x;
}

} // namespace detail

namespace {

// FastICA on the rows of X (n_features × n_samples): whitening
// W0 = P * diag(1 / sqrt(lambda)) of the top-k covariance eigenpairs, then the
// rotation. Returns the unmixing matrix W0 * W (n_features × k), Julia's
// ica_fit.W.
Eigen::MatrixXd fast_ica(const Eigen::MatrixXd& X, int k) {
    const Eigen::MatrixXd W0 = whitening_matrix(detail::ica_whitening_dense(X, k));
    Eigen::VectorXd mean_vec = X.rowwise().mean();
    Eigen::MatrixXd Xc = X.colwise() - mean_vec;
    Eigen::MatrixXd Z = W0.transpose() * Xc; // k × n
    return W0 * fast_ica_rotation(Z);
}

// The same with the truncated whitening of a sparse X: O(nnz * k) per Lanczos
// step instead of O(n^3). Equal to the dense result up to rounding and the
// eigenvector signs, which change the FastICA start point.
Eigen::MatrixXd fast_ica(const SpMat& X, int k) {
    const Eigen::MatrixXd W0 = whitening_matrix(detail::ica_whitening_truncated(X, k));

    // Z = W0' * (X - mean 1') without forming the dense centred X.
    const Eigen::VectorXd shift = W0.transpose() * row_means(X);
    Eigen::MatrixXd Z(k, X.cols());
    parallel_for(0, X.cols(), 64, [&](std::int64_t j) {
        Eigen::VectorXd col = -shift;
        for (SpMat::InnerIterator it(X, static_cast<Eigen::Index>(j)); it; ++it)
            col.noalias() += it.value() * W0.row(it.row()).transpose();
        Z.col(j) = col;
    });
    return W0 * fast_ica_rotation(Z);
}

} // namespace

// Port of Julia's DataFrame-based wrapper for cluster_molecules_on_mrf:
// FastICA of the gene co-occurrence matrix initialises the EM, with the hash
// initialisation as fallback if ICA throws.
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

    int n_genes = 0;
    for (int g : genes) if (g > n_genes) n_genes = g;

    // Gene co-occurrence matrix (n_genes × n_genes, 0-based genes) and its
    // FastICA unmixing matrix (n_genes × n_clusters). Large panels use a
    // sparse matrix (on whole-transcriptome panels ~0.1 % of the gene pairs
    // are neighbours) and the truncated whitening, as the dense one is cubic
    // in the number of genes.
    const bool truncated = n_genes > dense_whitening_max_genes;
    Eigen::MatrixXd cor_mat;
    SpMat sparse_cor_mat;
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
        Eigen::MatrixXd W = truncated ? fast_ica(sparse_cor_mat, n_clusters) : fast_ica(cor_mat, n_clusters);
        // Julia: (abs.(ica_fit.W) ./ sum(abs.(ica_fit.W), dims=1))'
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
        spdlog::warn("ICA did not converge ({}), falling back to hash initialization.", e.what());
    } catch (...) {
        spdlog::warn("ICA failed, falling back to hash initialization."); // GCOVR_EXCL_LINE: nothing throws non-std exceptions
    }

    // Core EM with ICA init (or nullptr → hash fallback)
    return cluster_molecules_on_mrf(
        genes, adj_list, confidence,
        n_clusters, tol, mrf_weight, max_iters, verbose,
        exprs_init_ptr.get()
    );
}

} // namespace baysor
