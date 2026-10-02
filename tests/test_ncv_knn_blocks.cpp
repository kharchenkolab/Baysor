// NCV neighbourhood k-NN blocks (src/processing/data_processing/neighborhood_composition.cpp):
// with large k the queries are processed in several byte-bounded blocks and,
// above k = 32, with the heap result set. The assembled count matrix and the
// projected vectors must be exactly what each query gives on its own.

#include <gtest/gtest.h>

#include "baysor/processing/data_processing/neighborhood_composition.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <random>
#include <vector>

namespace {

struct NcvInput {
    Eigen::MatrixXd pos;
    std::vector<int> genes;  // 1-based
    int n_genes = 0;
};

// Points on a coarse grid (many tied distances and exact duplicates) with
// 1-based genes.
NcvInput make_input(int n, int n_genes, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> coord(0, 60);
    std::uniform_int_distribution<int> gene(1, n_genes);
    NcvInput in;
    in.pos.resize(2, n);
    in.genes.resize(n);
    for (int i = 0; i < n; ++i) {
        in.pos(0, i) = 0.5 * coord(rng);
        in.pos(1, i) = 0.5 * coord(rng);
        in.genes[i] = gene(rng);
    }
    in.n_genes = n_genes;
    return in;
}

} // namespace

TEST(NcvKnnBlocks, CountMatrixOverSeveralBlocksMatchesPerQueryResults) {
    const NcvInput in = make_input(5000, 300, 5);
    const double floor = baysor::neighborhood_distance_floor(in.pos);
    std::vector<int> all_ids(in.pos.cols());
    for (int i = 0; i < in.pos.cols(); ++i) all_ids[i] = i;

    for (int k : {20, 600, 2000}) {  // one block / two blocks / three blocks
        Eigen::SparseMatrix<float> full = baysor::neighborhood_count_matrix_subset(
            in.pos, in.genes, all_ids, k, in.n_genes, nullptr, true, true, floor);
        ASSERT_EQ(full.cols(), in.pos.cols());
        ASSERT_TRUE(full.isCompressed());
        for (int q : {0, 1, 2047, 2048, 2049, 4095, 4096, 4659, 4660, 4999}) {
            Eigen::SparseMatrix<float> one = baysor::neighborhood_count_matrix_subset(
                in.pos, in.genes, std::vector<int>{q}, k, in.n_genes, nullptr, true, true, floor);
            std::vector<std::pair<int, float>> a, b;
            for (Eigen::SparseMatrix<float>::InnerIterator it(full, q); it; ++it) a.emplace_back(it.row(), it.value());
            for (Eigen::SparseMatrix<float>::InnerIterator it(one, 0); it; ++it) b.emplace_back(it.row(), it.value());
            ASSERT_FALSE(a.empty());
            ASSERT_EQ(a, b) << "k=" << k << " query=" << q;
        }
    }
}

TEST(NcvKnnBlocks, ProjectedVectorsOverSeveralBlocksMatchPerQueryResults) {
    const NcvInput in = make_input(5000, 300, 9);
    const double floor = baysor::neighborhood_distance_floor(in.pos);
    std::mt19937 rng(1);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    Eigen::MatrixXf gene_emb_t(8, in.n_genes);
    for (int g = 0; g < in.n_genes; ++g)
        for (int d = 0; d < 8; ++d) gene_emb_t(d, g) = nd(rng);

    for (int k : {600, 2000}) {
        Eigen::MatrixXf full = baysor::project_neighborhood_vectors(
            in.pos, in.genes, k, gene_emb_t, in.n_genes, nullptr, nullptr, true, true, floor, true);
        ASSERT_EQ(full.cols(), in.pos.cols());
        for (int q : {0, 2047, 2048, 4096, 4660, 4999}) {
            std::vector<int> ids{q};
            Eigen::MatrixXf one = baysor::project_neighborhood_vectors(
                in.pos, in.genes, k, gene_emb_t, in.n_genes, &ids, nullptr, true, true, floor, true);
            for (int d = 0; d < gene_emb_t.rows(); ++d)
                ASSERT_EQ(full(d, q), one(d, 0)) << "k=" << k << " query=" << q << " dim=" << d;
        }
    }
}
