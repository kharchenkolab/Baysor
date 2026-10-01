#include <gtest/gtest.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/file_reader.h>
#include <parquet/arrow/reader.h>
#include "baysor/data_loading/data.h"
#include "baysor/data_loading/prior_segmentation.h"
#include "baysor/utils/general.h"
#include "baysor/utils/options.h"
#include "baysor/utils/xenium.h"
#include "baysor/processing/utils/utils.h"
#include "baysor/processing/data_processing/triangulation.h"
#include "baysor/processing/data_processing/boundary_estimation.h"
#include "baysor/processing/data_processing/initialization.h"
#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/processing/data_processing/neighborhood_composition.h"
#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/models/bmm_data.h"
#include "baysor/processing/models/component.h"
#include "baysor/processing/distributions/mv_normal.h"
#include "baysor/processing/distributions/categorical_smoothed.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/reporting/output.h"
#include "baysor/reporting/run_report.h"

#include <Eigen/Dense>
#include "baysor/utils/thread_pool.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <parquet/arrow/writer.h>
#include <hdf5.h>
#include <tiffio.h>

// Helper: write a temp CSV file and return its path
static std::string write_temp_csv(const std::string& content, const std::string& suffix = ".csv") {
    char tmpl[] = "/tmp/baysor_test_XXXXXX";
    int fd = mkstemp(tmpl);
    EXPECT_GE(fd, 0);
    close(fd);
    std::string path = std::string(tmpl) + suffix;
    std::rename(tmpl, path.c_str());
    std::ofstream f(path);
    f << content;
    f.close();
    return path;
}

template <typename T>
static T arrow_test_unwrap(arrow::Result<T>&& result) {
    EXPECT_TRUE(result.ok()) << result.status().ToString();
    return std::move(*result);
}

static void arrow_expect_ok(const arrow::Status& status) {
    EXPECT_TRUE(status.ok()) << status.ToString();
}

static std::shared_ptr<arrow::Array> make_double_array(const std::vector<double>& values) {
    arrow::DoubleBuilder builder;
    arrow_expect_ok(builder.AppendValues(values));
    return arrow_test_unwrap(builder.Finish());
}

static std::shared_ptr<arrow::Array> make_string_array(const std::vector<std::string>& values) {
    arrow::StringBuilder builder;
    for (const auto& v : values) arrow_expect_ok(builder.Append(v));
    return arrow_test_unwrap(builder.Finish());
}

static std::shared_ptr<arrow::Array> make_int32_array(const std::vector<int32_t>& values) {
    arrow::Int32Builder builder;
    arrow_expect_ok(builder.AppendValues(values));
    return arrow_test_unwrap(builder.Finish());
}

static std::string write_temp_parquet(
    const std::vector<std::string>& names,
    const std::vector<std::shared_ptr<arrow::Array>>& arrays,
    int64_t row_group_size = 3
) {
    char tmpl[] = "/tmp/baysor_test_parquet_XXXXXX";
    int fd = mkstemp(tmpl);
    EXPECT_GE(fd, 0);
    close(fd);
    std::string path = std::string(tmpl) + ".parquet";
    std::rename(tmpl, path.c_str());

    std::vector<std::shared_ptr<arrow::Field>> fields;
    for (size_t i = 0; i < names.size(); ++i) {
        fields.push_back(arrow::field(names[i], arrays[i]->type()));
    }

    auto schema = arrow::schema(fields);
    auto table = arrow::Table::Make(schema, arrays);
    auto sink = arrow_test_unwrap(arrow::io::FileOutputStream::Open(path));
    arrow_expect_ok(parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), sink, row_group_size));
    arrow_expect_ok(sink->Close());
    return path;
}

static std::string write_temp_tiff_mask_u8(
    const std::vector<uint8_t>& pixels,
    uint32_t width,
    uint32_t height
) {
    char tmpl[] = "/tmp/baysor_test_mask_XXXXXX";
    int fd = mkstemp(tmpl);
    EXPECT_GE(fd, 0);
    close(fd);
    std::string path = std::string(tmpl) + ".tif";
    std::rename(tmpl, path.c_str());

    TIFF* tif = TIFFOpen(path.c_str(), "w");
    EXPECT_NE(tif, nullptr);
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, height);

    for (uint32_t row = 0; row < height; ++row) {
        auto* row_ptr = const_cast<uint8_t*>(pixels.data() + static_cast<size_t>(row) * width);
        EXPECT_EQ(TIFFWriteScanline(tif, row_ptr, row, 0), 1);
    }
    TIFFClose(tif);
    return path;
}

static baysor::BmmData<2> make_disconnected_bmm_data() {
    using baysor::AdjList;
    using baysor::BmmData;
    using baysor::CategoricalSmoothed;
    using baysor::Component;
    using baysor::MvNormal;
    using baysor::ShapePrior;

    BmmData<2> data;

    // Three disconnected spatial groups with sizes 3, 3, and 2.
    data.position_data.resize(2, 8);
    data.position_data <<
        0.0,  0.1,  0.2,   10.0, 10.1, 10.2,   20.0, 20.1,
        0.0,  0.0,  0.1,    0.0,  0.0,  0.1,    0.0,  0.0;

    // Internal representation is 0-based gene IDs.
    data.composition_data = {0, 0, 0, 1, 1, 1, 0, 0};
    data.confidence.assign(8, 1.0);  // disable the noise class in the E-step

    const int edge_src[] = {0, 1, 3, 4, 6};
    const int edge_dst[] = {1, 2, 4, 5, 7};
    const double edge_wt[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    data.adj_list = AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 5, 8);

    ShapePrior<2> prior;
    prior.std_values << 0.25, 0.25;
    prior.std_value_stds << 0.05, 0.05;
    prior.n_samples = 3;

    const Eigen::Matrix2d sigma = Eigen::Matrix2d::Identity() * 0.05;
    const Eigen::Vector2d centers[] = {
        (Eigen::Vector2d() << 0.1, 0.03).finished(),
        (Eigen::Vector2d() << 10.1, 0.03).finished(),
        (Eigen::Vector2d() << 20.05, 0.0).finished()
    };

    for (int ci = 0; ci < 3; ++ci) {
        MvNormal<2> pos_params(centers[ci], sigma);
        CategoricalSmoothed comp_params(2, 1.0);
        comp_params.set_uniform_counts(1.0f);
        data.components.emplace_back(pos_params, comp_params, prior, ci + 1);
    }

    data.assignment = {1, 1, 1, 2, 2, 2, 3, 3};
    data.max_component_guid = 3;
    data.noise_position_density = 1e-6;
    data.noise_density = 1e-6;
    data.prior_seg_confidence = 0.2;
    data.cluster_penalty_mult = 0.25;
    data.use_gene_smoothing = true;
    data.min_nuclei_frac = 0.1;
    data.mrf_strength = 0.1;
    data.real_edge_weight = 1.0;

    return data;
}

static baysor::AdjList make_chain_adj_list(int n_molecules) {
    std::vector<int> src;
    std::vector<int> dst;
    std::vector<double> wts;
    src.reserve(std::max(0, n_molecules - 1));
    dst.reserve(std::max(0, n_molecules - 1));
    wts.reserve(std::max(0, n_molecules - 1));

    for (int i = 0; i < n_molecules - 1; ++i) {
        src.push_back(i);
        dst.push_back(i + 1);
        wts.push_back(1.0 + 0.2 * (i % 3));
    }

    return baysor::AdjList::from_edge_list(
        src.data(), dst.data(), wts.data(),
        static_cast<int>(src.size()), n_molecules);
}

template <typename TActual, typename TExpected>
static void expect_vector_near(
    const std::vector<TActual>& actual,
    const std::vector<TExpected>& expected,
    double tol
) {
    ASSERT_EQ(actual.size(), expected.size());
    for (size_t i = 0; i < actual.size(); ++i) {
        EXPECT_NEAR(actual[i], expected[i], tol) << "index=" << i;
    }
}

static Eigen::SparseMatrix<float> make_random_count_matrix(
    int n_genes,
    int n_molecules,
    int nnz_per_col,
    int seed
) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> gene_dist(0, n_genes - 1);
    std::uniform_real_distribution<float> weight_dist(0.05f, 1.0f);

    std::vector<Eigen::Triplet<float>> trips;
    trips.reserve(static_cast<size_t>(n_molecules) * nnz_per_col);
    std::vector<int> ids;
    ids.reserve(nnz_per_col);

    for (int col = 0; col < n_molecules; ++col) {
        ids.clear();
        std::unordered_set<int> chosen;
        chosen.reserve(static_cast<size_t>(nnz_per_col) * 2);
        while (static_cast<int>(ids.size()) < nnz_per_col) {
            int id = gene_dist(rng);
            if (chosen.insert(id).second) ids.push_back(id);
        }
        for (int row : ids) {
            trips.emplace_back(row, col, std::log(weight_dist(rng) * 10000.0f + 1e-5f));
        }
    }

    Eigen::SparseMatrix<float> mat(n_genes, n_molecules);
    mat.setFromTriplets(trips.begin(), trips.end(),
                        [](float a, float b) { return a + b; });
    mat.makeCompressed();
    return mat;
}

static Eigen::MatrixXf estimate_gene_vectors_reference(
    const Eigen::SparseMatrix<float>& count_matrix,
    int n_components,
    bool per_molecule
) {
    int n_genes = static_cast<int>(count_matrix.rows());
    int n_mols = static_cast<int>(count_matrix.cols());
    if (n_genes == 0 || n_mols == 0) return Eigen::MatrixXf();

    std::mt19937 rng(42);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    Eigen::MatrixXf random_vectors(n_genes, n_components);
    for (int i = 0; i < n_genes; ++i)
        for (int j = 0; j < n_components; ++j)
            random_vectors(i, j) = normal(rng);

    Eigen::SparseMatrix<float> coexpr_sp = count_matrix * count_matrix.transpose();
    Eigen::MatrixXf coexpr(coexpr_sp);

    constexpr float var_clip = 0.05f;
    if (var_clip > 0 && n_genes > 1) {
        Eigen::VectorXf diag_vals = coexpr.diagonal();
        Eigen::VectorXf total_var = coexpr.rowwise().sum();

        std::vector<float> df_sorted(n_genes);
        for (int i = 0; i < n_genes; ++i)
            df_sorted[i] = (total_var(i) > 0) ? diag_vals(i) / total_var(i) : 0.0f;
        std::sort(df_sorted.begin(), df_sorted.end());
        float q = df_sorted[static_cast<int>((1.0f - var_clip) * (n_genes - 1))];

        for (int i = 0; i < n_genes; ++i)
            coexpr(i, i) = std::min(q * total_var(i), diag_vals(i));
    }

    Eigen::VectorXf row_sums = coexpr.rowwise().sum();
    Eigen::MatrixXf gene_emb = coexpr * random_vectors;
    for (int i = 0; i < n_genes; ++i) {
        if (row_sums(i) > 0) gene_emb.row(i) /= row_sums(i);
    }

    if (per_molecule) return (count_matrix.transpose() * gene_emb).transpose();
    return gene_emb.transpose();
}

static std::string format_int_vector(const std::vector<int>& values) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << values[i];
    }
    oss << "]";
    return oss.str();
}

static void expect_matrix_near(
    const Eigen::MatrixXd& actual,
    const Eigen::MatrixXd& expected,
    double tol
) {
    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());
    for (int r = 0; r < actual.rows(); ++r) {
        for (int c = 0; c < actual.cols(); ++c) {
            EXPECT_NEAR(actual(r, c), expected(r, c), tol)
                << "row=" << r << " col=" << c;
        }
    }
}

static std::vector<int> swap_binary_labels(const std::vector<int>& assignment) {
    std::vector<int> swapped = assignment;
    for (int& a : swapped) {
        if (a == 1) a = 2;
        else if (a == 2) a = 1;
    }
    return swapped;
}

static baysor::Component<2> make_one_gene_component(
    double center_x,
    double center_y,
    double variance,
    double prior_probability,
    int n_samples,
    int guid
) {
    Eigen::Vector2d center;
    center << center_x, center_y;
    Eigen::Matrix2d sigma = Eigen::Matrix2d::Identity() * variance;
    baysor::MvNormal<2> pos_params(center, sigma);
    baysor::CategoricalSmoothed comp_params(1, 1.0);
    comp_params.set_dense_counts({1.0f});
    baysor::Component<2> comp(pos_params, comp_params, std::nullopt, guid);
    comp.prior_probability = prior_probability;
    comp.n_samples = n_samples;
    comp.confidence = 1.0;
    return comp;
}

static baysor::Component<2> make_two_gene_component(
    double center_x,
    double center_y,
    double variance,
    double prior_probability,
    int n_samples,
    int guid
) {
    Eigen::Vector2d center;
    center << center_x, center_y;
    Eigen::Matrix2d sigma = Eigen::Matrix2d::Identity() * variance;
    baysor::MvNormal<2> pos_params(center, sigma);
    baysor::CategoricalSmoothed comp_params(2, 1.0);
    comp_params.set_dense_counts({1.0f, 1.0f});
    baysor::Component<2> comp(pos_params, comp_params, std::nullopt, guid);
    comp.prior_probability = prior_probability;
    comp.n_samples = n_samples;
    comp.confidence = 1.0;
    return comp;
}

static baysor::BmmData<2> make_two_component_competition_data() {
    baysor::BmmData<2> data;
    data.position_data.resize(2, 6);
    data.position_data <<
        0.15, 0.00, 0.00, 0.00, 0.15, 1.00,
        0.00, 0.00, 0.05, 0.10, 0.00, 0.00;
    data.composition_data = {0, 0, 0, 0, 0, 0};
    data.confidence.assign(6, 1.0);

    const int edge_src[] = {0, 0};
    const int edge_dst[] = {1, 4};
    const double edge_wt[] = {1.0, 1.0};
    data.adj_list = baysor::AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 2, 6);

    data.components.push_back(make_one_gene_component(
        /*center_x=*/0.0, /*center_y=*/0.0, /*variance=*/0.04,
        /*prior_probability=*/1.0, /*n_samples=*/3, /*guid=*/1));
    data.components.push_back(make_one_gene_component(
        /*center_x=*/0.15, /*center_y=*/0.0, /*variance=*/0.04,
        /*prior_probability=*/1.0, /*n_samples=*/1, /*guid=*/2));
    data.components.push_back(make_one_gene_component(
        /*center_x=*/1.0, /*center_y=*/0.0, /*variance=*/0.04,
        /*prior_probability=*/1.0, /*n_samples=*/1, /*guid=*/3));

    data.assignment = {0, 1, 1, 1, 2, 3};
    data.max_component_guid = 3;
    data.cluster_penalty_mult = 0.25;
    data.use_gene_smoothing = true;
    data.mrf_strength = 0.1;
    data.real_edge_weight = 1.0;
    return data;
}

static baysor::BmmData<2> make_last_component_competition_data() {
    baysor::BmmData<2> data;
    data.position_data.resize(2, 6);
    data.position_data <<
        0.15, 0.00, 0.00, 0.00, 0.15, 1.00,
        0.00, 0.00, 0.05, 0.10, 0.00, 0.00;
    data.composition_data = {0, 0, 0, 0, 0, 0};
    data.confidence.assign(6, 1.0);

    const int edge_src[] = {0, 0};
    const int edge_dst[] = {1, 4};
    const double edge_wt[] = {1.0, 1.0};
    data.adj_list = baysor::AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 2, 6);

    data.components.push_back(make_one_gene_component(
        /*center_x=*/0.0, /*center_y=*/0.0, /*variance=*/0.04,
        /*prior_probability=*/1.0, /*n_samples=*/3, /*guid=*/1));
    data.components.push_back(make_one_gene_component(
        /*center_x=*/1.0, /*center_y=*/0.0, /*variance=*/0.04,
        /*prior_probability=*/1.0, /*n_samples=*/1, /*guid=*/2));
    data.components.push_back(make_one_gene_component(
        /*center_x=*/0.15, /*center_y=*/0.0, /*variance=*/0.04,
        /*prior_probability=*/1.0, /*n_samples=*/1, /*guid=*/3));

    data.assignment = {0, 1, 1, 1, 3, 2};
    data.max_component_guid = 3;
    data.cluster_penalty_mult = 0.25;
    data.use_gene_smoothing = true;
    data.mrf_strength = 0.1;
    data.real_edge_weight = 1.0;
    return data;
}

static baysor::BmmData<2> make_noise_competition_data() {
    baysor::BmmData<2> data;
    data.position_data.resize(2, 2);
    data.position_data <<
        0.15, 0.15,
        0.00, 0.00;
    data.composition_data = {0, 0};
    data.confidence = {0.2, 1.0};

    const int edge_src[] = {0};
    const int edge_dst[] = {1};
    const double edge_wt[] = {1.0};
    data.adj_list = baysor::AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 1, 2);

    data.components.push_back(make_one_gene_component(
        /*center_x=*/0.15, /*center_y=*/0.0, /*variance=*/0.04,
        /*prior_probability=*/1e-3, /*n_samples=*/1, /*guid=*/1));

    data.assignment = {0, 1};
    data.max_component_guid = 1;
    data.use_gene_smoothing = true;
    data.mrf_strength = 0.1;
    data.real_edge_weight = 1.0;
    return data;
}

static baysor::BmmData<2> make_noise_density_skip_empty_component_data() {
    baysor::BmmData<2> data;
    data.position_data.resize(2, 2);
    data.position_data <<
        0.0, 0.1,
        0.0, 0.0;
    data.composition_data = {0, 1};
    data.confidence = {1.0, 1.0};

    data.components.push_back(make_two_gene_component(
        /*center_x=*/0.0, /*center_y=*/0.0, /*variance=*/0.04,
        /*prior_probability=*/1.0, /*n_samples=*/2, /*guid=*/1));
    data.components.push_back(make_two_gene_component(
        /*center_x=*/1.0, /*center_y=*/0.0, /*variance=*/0.04,
        /*prior_probability=*/1.0, /*n_samples=*/0, /*guid=*/2));

    data.assignment = {1, 1};
    data.max_component_guid = 2;
    data.noise_position_density = 2.0;
    return data;
}

// ============================================================================
// Gene encoding
// ============================================================================

TEST(GeneEncoding, BasicEncoding) {
    baysor::MoleculeData data;
    std::vector<std::string> genes = {"GeneC", "GeneA", "GeneB", "GeneA", "GeneC"};
    baysor::encode_genes(data, genes);

    // gene_names should be sorted alphabetically
    ASSERT_EQ(data.gene_names.size(), 3u);
    EXPECT_EQ(data.gene_names[0], "GeneA");
    EXPECT_EQ(data.gene_names[1], "GeneB");
    EXPECT_EQ(data.gene_names[2], "GeneC");

    // gene IDs should be 1-based
    ASSERT_EQ(data.gene.size(), 5u);
    EXPECT_EQ(data.gene[0], 3);  // GeneC -> 3
    EXPECT_EQ(data.gene[1], 1);  // GeneA -> 1
    EXPECT_EQ(data.gene[2], 2);  // GeneB -> 2
    EXPECT_EQ(data.gene[3], 1);  // GeneA -> 1
    EXPECT_EQ(data.gene[4], 3);  // GeneC -> 3
}

TEST(GeneEncoding, EmptyGenes) {
    baysor::MoleculeData data;
    std::vector<std::string> genes = {"A", "", "B"};
    baysor::encode_genes(data, genes);

    // Empty string is skipped from gene_names
    ASSERT_EQ(data.gene_names.size(), 2u);
    EXPECT_EQ(data.gene_names[0], "A");
    EXPECT_EQ(data.gene_names[1], "B");

    // Empty string maps to 0
    EXPECT_EQ(data.gene[0], 1);
    EXPECT_EQ(data.gene[1], 0);
    EXPECT_EQ(data.gene[2], 2);
}

// ============================================================================
// Gene filtering by count
// ============================================================================

TEST(GeneFiltering, FilterByCount) {
    baysor::MoleculeData data;
    data.x = {1, 2, 3, 4, 5, 6, 7};
    data.y = {1, 2, 3, 4, 5, 6, 7};
    std::vector<std::string> genes = {"A", "A", "A", "B", "C", "C", "C"};
    baysor::encode_genes(data, genes);

    // A: 3, B: 1, C: 3 — filter with min=2 should remove B
    baysor::filter_genes_by_count(data, 2);

    EXPECT_EQ(data.n_molecules(), 6);
    EXPECT_EQ(data.n_genes(), 2);
    EXPECT_EQ(data.gene_names[0], "A");
    EXPECT_EQ(data.gene_names[1], "C");
}

TEST(GeneFiltering, FilterByCountNoOp) {
    baysor::MoleculeData data;
    data.x = {1, 2, 3};
    data.y = {1, 2, 3};
    std::vector<std::string> genes = {"A", "B", "A"};
    baysor::encode_genes(data, genes);

    // A: 2, B: 1 — min=1 should not filter anything
    baysor::filter_genes_by_count(data, 1);
    EXPECT_EQ(data.n_molecules(), 3);
    EXPECT_EQ(data.n_genes(), 2);
}

TEST(GeneFiltering, FilterByCountCompactsMetadata) {
    baysor::MoleculeData data;
    data.x = {1, 2, 3, 4};
    data.y = {1, 2, 3, 4};
    data.confidence = {0.9, 0.8, 0.7, 0.6};
    std::vector<std::string> genes = {"A", "B", "A", "A"};
    baysor::encode_genes(data, genes);

    // A: 3, B: 1 — filter min=2 removes molecule 1 (index 1, gene B)
    baysor::filter_genes_by_count(data, 2);

    EXPECT_EQ(data.n_molecules(), 3);
    ASSERT_EQ(data.confidence.size(), 3u);
    EXPECT_DOUBLE_EQ(data.confidence[0], 0.9);
    EXPECT_DOUBLE_EQ(data.confidence[1], 0.7);
    EXPECT_DOUBLE_EQ(data.confidence[2], 0.6);
}

// ============================================================================
// Gene filtering by pattern
// ============================================================================

TEST(GeneFiltering, FilterByPattern) {
    baysor::MoleculeData data;
    data.x = {1, 2, 3, 4, 5};
    data.y = {1, 2, 3, 4, 5};
    std::vector<std::string> genes = {"Blank1", "GeneA", "Blank2", "GeneB", "MALAT1"};
    baysor::encode_genes(data, genes);

    // Filter genes matching "Blank*" pattern
    baysor::filter_genes_by_pattern(data, {"Blank*"});

    EXPECT_EQ(data.n_molecules(), 3);
    // Remaining genes: GeneA, GeneB, MALAT1
    EXPECT_EQ(data.n_genes(), 3);
}

TEST(GeneFiltering, FilterByMultiplePatterns) {
    baysor::MoleculeData data;
    data.x = {1, 2, 3, 4, 5};
    data.y = {1, 2, 3, 4, 5};
    std::vector<std::string> genes = {"Blank1", "GeneA", "Blank2", "GeneB", "MALAT1"};
    baysor::encode_genes(data, genes);

    baysor::filter_genes_by_pattern(data, {"Blank*", "MALAT1"});

    EXPECT_EQ(data.n_molecules(), 2);
    EXPECT_EQ(data.n_genes(), 2);
}

// ============================================================================
// Position matrix
// ============================================================================

TEST(MoleculeData, PositionMatrix2D) {
    baysor::MoleculeData data;
    data.x = {1.0, 2.0, 3.0};
    data.y = {4.0, 5.0, 6.0};

    auto mat = data.position_matrix();
    EXPECT_EQ(mat.rows(), 2);
    EXPECT_EQ(mat.cols(), 3);
    EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(mat(1, 2), 6.0);
}

TEST(MoleculeData, PositionMatrix3D) {
    baysor::MoleculeData data;
    data.x = {1.0, 2.0};
    data.y = {3.0, 4.0};
    data.z = {5.0, 6.0};

    auto mat = data.position_matrix();
    EXPECT_EQ(mat.rows(), 3);
    EXPECT_EQ(mat.cols(), 2);
    EXPECT_DOUBLE_EQ(mat(2, 1), 6.0);
}

TEST(MoleculeData, Dimensionality) {
    baysor::MoleculeData data;
    data.x = {1.0};
    data.y = {2.0};
    EXPECT_FALSE(data.is_3d());
    EXPECT_EQ(data.n_dims(), 2);

    data.z = {3.0};
    EXPECT_TRUE(data.is_3d());
    EXPECT_EQ(data.n_dims(), 3);
}

TEST(MvNormal, JuliaParity2DNormalization) {
    Eigen::Vector2d mu;
    mu << 0.0, 0.0;
    Eigen::Matrix2d sigma = Eigen::Matrix2d::Identity();
    baysor::MvNormal<2> dist(mu, sigma);

    const double x[2] = {0.0, 0.0};
    EXPECT_NEAR(dist.log_pdf(x), -2.756815599614018, 1e-12);
    EXPECT_NEAR(dist.pdf(x), 0.06349363593424098, 1e-12);
}

TEST(MvNormal, JuliaParityShapePriorSignedOffDiagonalZeroing) {
    Eigen::Matrix2d sigma;
    sigma << 4.0, -2.0,
            -2.0, 3.0;

    baysor::ShapePrior<2> prior;
    prior.std_values << 2.0, 2.0;
    prior.std_value_stds << 0.5, 0.5;
    prior.n_samples = 20;

    baysor::adjust_cov_by_prior<2>(sigma, prior, /*n_samples=*/10);

    // Julia zeros the signed-negative off-diagonal entry before the
    // eigendecomposition, so the posterior covariance stays diagonal here.
    EXPECT_NEAR(sigma(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(sigma(1, 0), 0.0, 1e-12);
}

TEST(Component, IndexedMaximizeMatchesContiguousPath) {
    Eigen::Vector2d mu;
    mu << 0.0, 0.0;
    Eigen::Matrix2d sigma = Eigen::Matrix2d::Identity();

    baysor::ShapePrior<2> prior;
    prior.std_values << 1.5, 1.5;
    prior.std_value_stds << 0.2, 0.2;
    prior.n_samples = 5;

    baysor::CategoricalSmoothed comp_params(2, 1.0);
    baysor::Component<2> contiguous(
        baysor::MvNormal<2>(mu, sigma), comp_params, prior, /*guid=*/1);
    baysor::Component<2> indexed(
        baysor::MvNormal<2>(mu, sigma), comp_params, prior, /*guid=*/1);

    Eigen::MatrixXd pos_data(2, 4);
    pos_data <<
        0.0, 2.0, 4.0, 5.0,
        0.0, 1.0, 1.5, 3.0;
    std::vector<int> gene_ids = {0, 1, 1, 0};
    std::vector<double> nuclei_probs = {0.1, 0.2, 0.8, 0.9};
    const int mol_ids[] = {0, 2, 3};

    double pos_buf[] = {
        0.0, 0.0,
        4.0, 1.5,
        5.0, 3.0
    };
    int gene_buf[] = {0, 1, 0};
    double nuc_buf[] = {0.1, 0.8, 0.9};

    contiguous.maximize(pos_buf, /*stride=*/2, gene_buf, /*n_points=*/3, nuc_buf,
                        /*min_nuclei_frac=*/0.2, /*freeze_position=*/false,
                        /*freeze_composition=*/false);
    indexed.maximize_indexed(pos_data, gene_ids, mol_ids, /*n_points=*/3, &nuclei_probs,
                             /*min_nuclei_frac=*/0.2, /*freeze_position=*/false,
                             /*freeze_composition=*/false);

    EXPECT_EQ(indexed.n_samples, contiguous.n_samples);
    EXPECT_EQ(indexed.composition_params.n_genes, contiguous.composition_params.n_genes);
    EXPECT_NEAR(indexed.composition_params.sum_counts, contiguous.composition_params.sum_counts, 1e-12);
    expect_vector_near(indexed.composition_params.dense_counts(), contiguous.composition_params.dense_counts(), 1e-6);
    EXPECT_NEAR(indexed.position_params.mu(0), contiguous.position_params.mu(0), 1e-12);
    EXPECT_NEAR(indexed.position_params.mu(1), contiguous.position_params.mu(1), 1e-12);
    EXPECT_NEAR(indexed.position_params.sigma(0, 0), contiguous.position_params.sigma(0, 0), 1e-12);
    EXPECT_NEAR(indexed.position_params.sigma(0, 1), contiguous.position_params.sigma(0, 1), 1e-12);
    EXPECT_NEAR(indexed.position_params.sigma(1, 0), contiguous.position_params.sigma(1, 0), 1e-12);
    EXPECT_NEAR(indexed.position_params.sigma(1, 1), contiguous.position_params.sigma(1, 1), 1e-12);
    EXPECT_NEAR(indexed.confidence, contiguous.confidence, 1e-12);
    EXPECT_NEAR(indexed.confidence, 0.86, 1e-12);
}

// ============================================================================
// CSV loading
// ============================================================================

TEST(DataLoading, LoadCSVBasic) {
    auto path = write_temp_csv(
        "x,y,gene\n"
        "1.0,2.0,A\n"
        "3.0,4.0,B\n"
        "5.0,6.0,A\n"
    );

    baysor::MoleculeInputOptions opts;
    auto data = baysor::load_molecules(path, opts);

    EXPECT_EQ(data.n_molecules(), 3);
    EXPECT_EQ(data.n_genes(), 2);
    EXPECT_DOUBLE_EQ(data.x[0], 1.0);
    EXPECT_DOUBLE_EQ(data.y[2], 6.0);
    EXPECT_FALSE(data.is_3d());

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVWith3D) {
    auto path = write_temp_csv(
        "x,y,z,gene\n"
        "1.0,2.0,10.0,A\n"
        "3.0,4.0,20.0,B\n"
    );

    baysor::MoleculeInputOptions opts;
    auto data = baysor::load_molecules(path, opts);

    EXPECT_TRUE(data.is_3d());
    EXPECT_EQ(data.n_dims(), 3);
    EXPECT_DOUBLE_EQ(data.z[1], 20.0);

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVDropConstantZ) {
    auto path = write_temp_csv(
        "x,y,z,gene\n"
        "1.0,2.0,0.0,A\n"
        "3.0,4.0,0.0,B\n"
    );

    baysor::MoleculeInputOptions opts;
    auto data = baysor::load_molecules(path, opts);

    // z is constant -> should be dropped
    EXPECT_FALSE(data.is_3d());
    EXPECT_TRUE(data.z.empty());

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVForce2D) {
    auto path = write_temp_csv(
        "x,y,z,gene\n"
        "1.0,2.0,10.0,A\n"
        "3.0,4.0,20.0,B\n"
    );

    baysor::MoleculeInputOptions opts;
    opts.force_2d = true;
    auto data = baysor::load_molecules(path, opts);

    EXPECT_FALSE(data.is_3d());
    EXPECT_TRUE(data.z.empty());

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVWithConfidence) {
    auto path = write_temp_csv(
        "x,y,gene,confidence\n"
        "1.0,2.0,A,0.95\n"
        "3.0,4.0,B,0.80\n"
    );

    baysor::MoleculeInputOptions opts;
    auto data = baysor::load_molecules(path, opts);

    ASSERT_EQ(data.confidence.size(), 2u);
    EXPECT_DOUBLE_EQ(data.confidence[0], 0.95);
    EXPECT_DOUBLE_EQ(data.confidence[1], 0.80);

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVCustomColumns) {
    auto path = write_temp_csv(
        "pos_x,pos_y,transcript\n"
        "1.0,2.0,GeneA\n"
        "3.0,4.0,GeneB\n"
    );

    baysor::MoleculeInputOptions opts;
    opts.x_col = "pos_x";
    opts.y_col = "pos_y";
    opts.gene_col = "transcript";
    auto data = baysor::load_molecules(path, opts);

    EXPECT_EQ(data.n_molecules(), 2);
    EXPECT_DOUBLE_EQ(data.x[0], 1.0);
    EXPECT_EQ(data.gene_names[0], "GeneA");

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVWithGeneFilter) {
    auto path = write_temp_csv(
        "x,y,gene\n"
        "1.0,1.0,A\n"
        "2.0,2.0,A\n"
        "3.0,3.0,A\n"
        "4.0,4.0,B\n"
        "5.0,5.0,C\n"
        "6.0,6.0,C\n"
        "7.0,7.0,C\n"
    );

    baysor::MoleculeInputOptions opts;
    opts.min_molecules_per_gene = 3;
    auto data = baysor::load_molecules(path, opts);

    // B has 1, C has 3, A has 3 — B filtered out
    EXPECT_EQ(data.n_molecules(), 6);
    EXPECT_EQ(data.n_genes(), 2);

    std::remove(path.c_str());
}

TEST(DataLoading, LoadCSVWithExcludeGenes) {
    auto path = write_temp_csv(
        "x,y,gene\n"
        "1.0,1.0,Blank1\n"
        "2.0,2.0,GeneA\n"
        "3.0,3.0,Blank2\n"
        "4.0,4.0,GeneB\n"
    );

    baysor::MoleculeInputOptions opts;
    opts.exclude_genes = "Blank*";
    auto data = baysor::load_molecules(path, opts);

    EXPECT_EQ(data.n_molecules(), 2);
    EXPECT_EQ(data.n_genes(), 2);

    std::remove(path.c_str());
}

TEST(DataLoading, LoadParquetWithTwoPassFiltersAndPriorColumn) {
    auto path = write_temp_parquet(
        {"x_location", "y_location", "z_location", "feature_name", "cell_id", "qv"},
        {
            make_double_array({0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 100.0, 6.0}),
            make_double_array({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}),
            make_double_array({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}),
            make_string_array({"GeneA", "GeneA", "GeneA", "GeneB", "GeneB", "NegControl1", "GeneA", "GeneC"}),
            make_string_array({"cell1", "cell1", "cell1", "cell2", "UNASSIGNED", "cell3", "cell1", "cell4"}),
            make_int32_array({30, 10, 30, 30, 30, 30, 30, 30})
        },
        /*row_group_size=*/3
    );

    baysor::MoleculeInputOptions opts;
    opts.x_col = "x_location";
    opts.y_col = "y_location";
    opts.z_col = "z_location";
    opts.gene_col = "feature_name";
    opts.qv_col = "qv";
    opts.min_qv = 20.0;
    opts.x_max = 10.0;
    opts.min_molecules_per_gene = 2;
    opts.exclude_genes = "NegControl*";

    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Column;
    prior.column_name = "cell_id";
    prior.unassigned_label = "UNASSIGNED";
    prior.min_molecules_per_segment = 2;

    auto data = baysor::load_molecules(path, opts, prior);

    EXPECT_EQ(data.n_molecules(), 4);
    EXPECT_EQ(data.gene_names.size(), 2u);
    EXPECT_EQ(data.gene_names[0], "GeneA");
    EXPECT_EQ(data.gene_names[1], "GeneB");
    EXPECT_EQ(data.gene, (std::vector<int>{1, 1, 2, 2}));
    EXPECT_TRUE(data.z.empty());
    EXPECT_EQ(data.prior_segmentation, (std::vector<int>{1, 1, 0, 0}));

    std::remove(path.c_str());
}

TEST(PriorSegmentation, EncodePriorLabels) {
    std::vector<std::string> labels = {"cell2", "cell1", "cell2", "UNASSIGNED", "cell1"};
    auto encoded = baysor::encode_prior_labels(labels, "UNASSIGNED", 0);
    EXPECT_EQ(encoded, (std::vector<int>{2, 1, 2, 0, 1}));
}

// ============================================================================
// Prior segmentation
// ============================================================================

TEST(PriorSegmentation, DetectType) {
    EXPECT_EQ(baysor::detect_prior_seg_type(""), baysor::PriorInputType::None);
    EXPECT_EQ(baysor::detect_prior_seg_type(":cell_id"), baysor::PriorInputType::Column);
    EXPECT_EQ(baysor::detect_prior_seg_type("/path/to/mask.tiff"), baysor::PriorInputType::Image);
    EXPECT_EQ(baysor::detect_prior_seg_type("/path/to/boundaries.csv"), baysor::PriorInputType::Boundary);
    EXPECT_EQ(baysor::detect_prior_seg_type("/path/to/boundaries.parquet"), baysor::PriorInputType::Boundary);
}

TEST(PriorSegmentation, ApplyInputSpecPreservesInterpretationOptions) {
    baysor::PriorInputOptions prior;
    prior.unassigned_label = "UNASSIGNED";
    prior.min_molecules_per_segment = 12;
    prior.estimate_scale_from_prior = false;

    baysor::apply_prior_input_spec(prior, ":cell_id");

    EXPECT_EQ(prior.type, baysor::PriorInputType::Column);
    EXPECT_EQ(prior.column_name, "cell_id");
    EXPECT_TRUE(prior.path.empty());
    EXPECT_EQ(prior.unassigned_label, "UNASSIGNED");
    EXPECT_EQ(prior.min_molecules_per_segment, 12);
    EXPECT_FALSE(prior.estimate_scale_from_prior);

    baysor::apply_prior_input_spec(prior, "/path/to/boundaries.parquet");

    EXPECT_EQ(prior.type, baysor::PriorInputType::Boundary);
    EXPECT_EQ(prior.path, "/path/to/boundaries.parquet");
    EXPECT_TRUE(prior.column_name.empty());
    EXPECT_EQ(prior.unassigned_label, "UNASSIGNED");
    EXPECT_EQ(prior.min_molecules_per_segment, 12);
    EXPECT_FALSE(prior.estimate_scale_from_prior);
}

TEST(PriorSegmentation, FilterLabels) {
    std::vector<int> labels = {1, 1, 1, 2, 3, 3, 3, 3};
    baysor::filter_segmentation_labels(labels, 3);

    // Segment 1: 3 molecules (passes), segment 2: 1 (filtered), segment 3: 4 (passes)
    EXPECT_EQ(labels[0], 1);
    EXPECT_EQ(labels[3], 0);  // was 2, now filtered
    EXPECT_EQ(labels[4], 3);
}

TEST(PriorSegmentation, FilterLabelsNoOp) {
    std::vector<int> labels = {1, 1, 2, 2};
    baysor::filter_segmentation_labels(labels, 2);

    // Both segments have 2 molecules, should pass
    EXPECT_EQ(labels[0], 1);
    EXPECT_EQ(labels[2], 2);
}

TEST(PriorSegmentation, ParseFromColumn) {
    auto path = write_temp_csv(
        "x,y,gene,cell_id\n"
        "1.0,1.0,A,cell1\n"
        "2.0,2.0,A,cell1\n"
        "3.0,3.0,A,cell1\n"
        "4.0,4.0,B,cell2\n"
        "5.0,5.0,B,0\n"
        "6.0,6.0,B,cell2\n"
    );

    auto labels = baysor::parse_prior_from_column(path, "cell_id", "0", 0);

    ASSERT_EQ(labels.size(), 6u);
    // cell1 -> 1, cell2 -> 2 (alphabetically sorted), "0" -> 0 (unassigned)
    EXPECT_EQ(labels[0], 1);  // cell1
    EXPECT_EQ(labels[1], 1);  // cell1
    EXPECT_EQ(labels[2], 1);  // cell1
    EXPECT_EQ(labels[3], 2);  // cell2
    EXPECT_EQ(labels[4], 0);  // unassigned
    EXPECT_EQ(labels[5], 2);  // cell2

    std::remove(path.c_str());
}

TEST(PriorSegmentation, ParseFromColumnWithFilter) {
    auto path = write_temp_csv(
        "x,y,gene,cell_id\n"
        "1.0,1.0,A,cell1\n"
        "2.0,2.0,A,cell1\n"
        "3.0,3.0,A,cell1\n"
        "4.0,4.0,B,cell2\n"
        "5.0,5.0,B,0\n"
    );

    // cell2 has only 1 molecule, filter with min=2
    auto labels = baysor::parse_prior_from_column(path, "cell_id", "0", 2);

    EXPECT_EQ(labels[3], 0);  // cell2 filtered out

    std::remove(path.c_str());
}

TEST(PriorSegmentation, EstimateScale) {
    // Create a simple grid of 4 cell centers spaced 10 units apart
    // Cells at (0,0), (10,0), (0,10), (10,10)
    // Each cell has 5 molecules clustered around its center
    Eigen::MatrixXd pos(2, 20);
    std::vector<int> assignment(20);

    // Cell 1: around (0, 0)
    for (int i = 0; i < 5; ++i) {
        pos(0, i) = 0.0 + i * 0.1;
        pos(1, i) = 0.0 + i * 0.1;
        assignment[i] = 1;
    }
    // Cell 2: around (10, 0)
    for (int i = 0; i < 5; ++i) {
        pos(0, 5 + i) = 10.0 + i * 0.1;
        pos(1, 5 + i) = 0.0 + i * 0.1;
        assignment[5 + i] = 2;
    }
    // Cell 3: around (0, 10)
    for (int i = 0; i < 5; ++i) {
        pos(0, 10 + i) = 0.0 + i * 0.1;
        pos(1, 10 + i) = 10.0 + i * 0.1;
        assignment[10 + i] = 3;
    }
    // Cell 4: around (10, 10)
    for (int i = 0; i < 5; ++i) {
        pos(0, 15 + i) = 10.0 + i * 0.1;
        pos(1, 15 + i) = 10.0 + i * 0.1;
        assignment[15 + i] = 4;
    }

    auto [scale, scale_std] = baysor::estimate_scale_from_assignment(pos, assignment, 3);

    // NN distance between centers ~10.0, so radius ~5.0
    EXPECT_NEAR(scale, 5.0, 0.5);
    EXPECT_GE(scale_std, 0.0);
}

TEST(PriorSegmentation, EstimateScaleMatchesExactNearestCenterReference) {
    constexpr int n_cells = 6;
    constexpr int mols_per_cell = 4;
    Eigen::MatrixXd centers(3, n_cells);
    centers.col(0) << 0.0, 0.0, 0.0;
    centers.col(1) << 5.0, 0.0, 0.0;
    centers.col(2) << 0.0, 7.0, 1.0;
    centers.col(3) << 10.0, 1.0, 2.0;
    centers.col(4) << 13.0, 5.0, 1.0;
    centers.col(5) << 20.0, 0.0, 3.0;

    Eigen::MatrixXd pos(3, n_cells * mols_per_cell);
    std::vector<int> assignment(n_cells * mols_per_cell);
    const Eigen::Vector3d offsets[mols_per_cell] = {
        Eigen::Vector3d(-0.1, -0.1, 0.0),
        Eigen::Vector3d( 0.1, -0.1, 0.0),
        Eigen::Vector3d(-0.1,  0.1, 0.0),
        Eigen::Vector3d( 0.1,  0.1, 0.0)
    };
    for (int c = 0; c < n_cells; ++c) {
        for (int m = 0; m < mols_per_cell; ++m) {
            const int idx = c * mols_per_cell + m;
            pos.col(idx) = centers.col(c) + offsets[m];
            assignment[idx] = c + 1;
        }
    }

    std::vector<double> radii(n_cells);
    for (int i = 0; i < n_cells; ++i) {
        double min_dist = std::numeric_limits<double>::max();
        for (int j = 0; j < n_cells; ++j) {
            if (i == j) continue;
            min_dist = std::min(min_dist, (centers.col(i) - centers.col(j)).norm());
        }
        radii[i] = min_dist / 2.0;
    }
    std::sort(radii.begin(), radii.end());
    const double expected_scale = (radii[n_cells / 2 - 1] + radii[n_cells / 2]) / 2.0;

    std::vector<double> abs_devs(n_cells);
    for (int i = 0; i < n_cells; ++i) {
        abs_devs[i] = std::abs(radii[i] - expected_scale);
    }
    std::sort(abs_devs.begin(), abs_devs.end());
    const double expected_std =
        ((abs_devs[n_cells / 2 - 1] + abs_devs[n_cells / 2]) / 2.0) * 1.4826;

    int old_threads = baysor::thread_pool_size();
    baysor::set_thread_pool_size(1);
    auto [scale_1, scale_std_1] = baysor::estimate_scale_from_assignment(pos, assignment, mols_per_cell);
    baysor::set_thread_pool_size(4);
    auto [scale_4, scale_std_4] = baysor::estimate_scale_from_assignment(pos, assignment, mols_per_cell);
    baysor::set_thread_pool_size(old_threads);

    EXPECT_NEAR(scale_1, expected_scale, 1e-10);
    EXPECT_NEAR(scale_std_1, expected_std, 1e-10);
    EXPECT_NEAR(scale_4, expected_scale, 1e-10);
    EXPECT_NEAR(scale_std_4, expected_std, 1e-10);
}

// ============================================================================
// Top-level prior segmentation loading
// ============================================================================

TEST(PriorSegmentation, LoadFromColumn) {
    auto path = write_temp_csv(
        "x,y,gene,cell_id\n"
        "1.0,1.0,A,cell1\n"
        "2.0,2.0,A,cell1\n"
        "3.0,3.0,A,cell1\n"
        "10.0,10.0,B,cell2\n"
        "11.0,11.0,B,cell2\n"
        "12.0,12.0,B,cell2\n"
        "20.0,20.0,A,cell3\n"
        "21.0,21.0,A,cell3\n"
        "22.0,22.0,A,cell3\n"
    );

    baysor::MoleculeInputOptions opts;
    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Column;
    prior.column_name = "cell_id";

    auto data = baysor::load_molecules(path, opts, prior);

    auto [scale, scale_std] = baysor::load_prior_segmentation(
        data, prior, 3);

    ASSERT_EQ(data.prior_segmentation.size(), 9u);
    EXPECT_GT(data.prior_segmentation[0], 0);  // assigned
    EXPECT_GT(scale, 0.0);

    std::remove(path.c_str());
}

TEST(PriorSegmentation, LoadFromBoundaryCsv) {
    auto coord_path = write_temp_csv(
        "x,y,gene\n"
        "1.0,1.0,A\n"
        "2.0,2.0,A\n"
        "4.0,4.0,A\n"
        "11.0,1.0,B\n"
        "12.0,2.0,B\n"
        "14.0,4.0,B\n"
        "8.0,2.0,C\n"
    );

    auto boundary_path = write_temp_csv(
        "cell_id,vertex_x,vertex_y\n"
        "cell1,0,0\n"
        "cell1,0,5\n"
        "cell1,5,5\n"
        "cell1,5,0\n"
        "cell2,10,0\n"
        "cell2,10,5\n"
        "cell2,15,5\n"
        "cell2,15,0\n"
    );

    baysor::MoleculeInputOptions opts;
    auto data = baysor::load_molecules(coord_path, opts);

    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Boundary;
    prior.path = boundary_path;

    auto [scale, scale_std] = baysor::load_prior_segmentation(
        data, prior, 2);

    EXPECT_EQ(data.prior_segmentation, (std::vector<int>{1, 1, 1, 2, 2, 2, 0}));
    EXPECT_LT(scale, 0.0);
    EXPECT_LT(scale_std, 0.0);

    std::remove(coord_path.c_str());
    std::remove(boundary_path.c_str());
}

TEST(PriorSegmentation, LoadNone) {
    baysor::MoleculeData data;
    data.x = {1.0};
    data.y = {1.0};

    auto [scale, scale_std] = baysor::load_prior_segmentation(
        data, baysor::PriorInputOptions{}, 3);

    EXPECT_DOUBLE_EQ(scale, -1.0);
    EXPECT_TRUE(data.prior_segmentation.empty());
}

TEST(PriorSegmentation, BinaryMaskScaleUsesFilteredComponentsOnly) {
    const uint32_t width = 20, height = 20;
    std::vector<uint8_t> pixels(width * height, 0);
    auto set_rect = [&](int r0, int r1, int c0, int c1) {
        for (int r = r0; r <= r1; ++r) {
            for (int c = c0; c <= c1; ++c) {
                pixels[static_cast<size_t>(r) * width + c] = 1;
            }
        }
    };

    // Three kept 2x2 components (area 4 each).
    set_rect(0, 1, 0, 1);
    set_rect(0, 1, 5, 6);
    set_rect(0, 1, 10, 11);
    // Three filtered 6x6 components (area 36 each).
    set_rect(9, 14, 0, 5);
    set_rect(9, 14, 7, 12);
    set_rect(9, 14, 14, 19);

    auto mask_path = write_temp_tiff_mask_u8(pixels, width, height);
    auto csv_path = write_temp_csv(
        "x,y,gene\n"
        "1,1,A\n"
        "2,2,A\n"
        "6,1,A\n"
        "7,2,A\n"
        "11,1,A\n"
        "12,2,A\n"
        "1,10,A\n"
        "8,10,A\n"
        "15,10,A\n"
    );

    baysor::MoleculeInputOptions opts;
    opts.min_molecules_per_cell = 2;
    baysor::fill_and_check_molecule_input_options(opts);
    auto data = baysor::load_molecules(csv_path, opts);

    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Image;
    prior.path = mask_path;
    prior.min_molecules_per_segment = 2;
    prior.estimate_scale_from_prior = true;

    auto [scale, scale_std] = baysor::load_prior_segmentation(
        data, prior, /*min_molecules_per_cell=*/2);

    const double expected_radius = std::sqrt(4.0 / baysor::kPi);
    EXPECT_NEAR(scale, expected_radius, 1e-6);
    EXPECT_NEAR(scale_std, 0.0, 1e-6);

    std::remove(mask_path.c_str());
    std::remove(csv_path.c_str());
}

TEST(PriorSegmentation, LabeledMaskPreservesTouchingLabelsAndUsesLabelAreas) {
    const uint32_t width = 9, height = 2;
    std::vector<uint8_t> pixels(width * height, 0);
    for (uint32_t row = 0; row < height; ++row) {
        for (uint32_t col = 0; col < 3; ++col) {
            pixels[static_cast<size_t>(row) * width + col] = 1;
            pixels[static_cast<size_t>(row) * width + col + 3] = 2;
            pixels[static_cast<size_t>(row) * width + col + 6] = 3;
        }
    }

    auto mask_path = write_temp_tiff_mask_u8(pixels, width, height);
    auto csv_path = write_temp_csv(
        "x,y,gene\n"
        "1,1,A\n"
        "2,2,A\n"
        "4,1,A\n"
        "5,2,A\n"
        "7,1,A\n"
        "8,2,A\n"
    );

    baysor::MoleculeInputOptions opts;
    opts.min_molecules_per_cell = 2;
    baysor::fill_and_check_molecule_input_options(opts);
    auto data = baysor::load_molecules(csv_path, opts);

    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Image;
    prior.path = mask_path;
    prior.min_molecules_per_segment = 2;
    prior.estimate_scale_from_prior = true;

    auto [scale, scale_std] = baysor::load_prior_segmentation(
        data, prior, /*min_molecules_per_cell=*/2);

    const std::vector<int> expected_labels = {1, 1, 2, 2, 3, 3};
    EXPECT_EQ(data.prior_segmentation, expected_labels);

    const double expected_radius = std::sqrt(6.0 / baysor::kPi);
    EXPECT_NEAR(scale, expected_radius, 1e-6);
    EXPECT_NEAR(scale_std, 0.0, 1e-6);

    std::remove(mask_path.c_str());
    std::remove(csv_path.c_str());
}

TEST(PriorSegmentation, BinaryMaskCropWindowPreservesAssignmentsAndAreas) {
    const uint32_t width = 120, height = 120;
    std::vector<uint8_t> pixels(width * height, 0);

    for (uint32_t r = 10; r <= 14; ++r) {
        for (uint32_t c = 20; c <= 24; ++c) {
            pixels[static_cast<size_t>(r) * width + c] = 1;
        }
    }
    for (uint32_t r = 20; r <= 24; ++r) {
        for (uint32_t c = 40; c <= 44; ++c) {
            pixels[static_cast<size_t>(r) * width + c] = 1;
        }
    }
    for (uint32_t r = 30; r <= 34; ++r) {
        for (uint32_t c = 60; c <= 64; ++c) {
            pixels[static_cast<size_t>(r) * width + c] = 1;
        }
    }

    auto mask_path = write_temp_tiff_mask_u8(pixels, width, height);
    auto csv_path = write_temp_csv(
        "x,y,gene\n"
        "21,11,A\n"
        "25,11,A\n"
        "21,15,A\n"
        "25,15,A\n"
        "41,21,A\n"
        "45,21,A\n"
        "41,25,A\n"
        "45,25,A\n"
        "61,31,A\n"
        "65,31,A\n"
        "61,35,A\n"
        "65,35,A\n"
    );

    baysor::MoleculeInputOptions opts;
    opts.min_molecules_per_cell = 2;
    baysor::fill_and_check_molecule_input_options(opts);
    auto data = baysor::load_molecules(csv_path, opts);

    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Image;
    prior.path = mask_path;
    prior.min_molecules_per_segment = 2;
    prior.estimate_scale_from_prior = true;

    auto [scale, scale_std] = baysor::load_prior_segmentation(data, prior, 2);

    EXPECT_EQ(data.prior_segmentation,
              (std::vector<int>{1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3}));
    EXPECT_NEAR(scale, std::sqrt(25.0 / baysor::kPi), 1e-6);
    EXPECT_NEAR(scale_std, 0.0, 1e-6);

    std::remove(mask_path.c_str());
    std::remove(csv_path.c_str());
}

TEST(PriorSegmentation, LabeledMaskCropWindowPreservesAssignmentsAndAreas) {
    const uint32_t width = 140, height = 140;
    std::vector<uint8_t> pixels(width * height, 0);

    for (uint32_t r = 20; r <= 22; ++r) {
        for (uint32_t c = 30; c <= 32; ++c) {
            pixels[static_cast<size_t>(r) * width + c] = 1;
        }
        for (uint32_t c = 35; c <= 37; ++c) {
            pixels[static_cast<size_t>(r) * width + c] = 2;
        }
        for (uint32_t c = 40; c <= 42; ++c) {
            pixels[static_cast<size_t>(r) * width + c] = 3;
        }
    }

    auto mask_path = write_temp_tiff_mask_u8(pixels, width, height);
    auto csv_path = write_temp_csv(
        "x,y,gene\n"
        "31,21,A\n"
        "33,23,A\n"
        "36,21,A\n"
        "38,23,A\n"
        "41,21,A\n"
        "43,23,A\n"
    );

    baysor::MoleculeInputOptions opts;
    opts.min_molecules_per_cell = 2;
    baysor::fill_and_check_molecule_input_options(opts);
    auto data = baysor::load_molecules(csv_path, opts);

    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Image;
    prior.path = mask_path;
    prior.min_molecules_per_segment = 2;
    prior.estimate_scale_from_prior = true;

    auto [scale, scale_std] = baysor::load_prior_segmentation(data, prior, 2);

    EXPECT_EQ(data.prior_segmentation, (std::vector<int>{1, 1, 2, 2, 3, 3}));
    EXPECT_NEAR(scale, std::sqrt(9.0 / baysor::kPi), 1e-6);
    EXPECT_NEAR(scale_std, 0.0, 1e-6);

    std::remove(mask_path.c_str());
    std::remove(csv_path.c_str());
}

// ============================================================================
// Options
// ============================================================================

TEST(Options, ParseScaleStd) {
    EXPECT_DOUBLE_EQ(baysor::parse_scale_std("25%", 10.0), 2.5);
    EXPECT_DOUBLE_EQ(baysor::parse_scale_std("50%", 10.0), 5.0);
    EXPECT_DOUBLE_EQ(baysor::parse_scale_std("3.0", 10.0), 3.0);
    EXPECT_DOUBLE_EQ(baysor::parse_scale_std("", 10.0), 2.5);  // default 25%
}

TEST(Options, DefaultParamValues) {
    EXPECT_EQ(baysor::default_param_value("min_molecules_per_segment", 12), 3);
    EXPECT_EQ(baysor::default_param_value("confidence_nn_id", 10), 6);
    EXPECT_EQ(baysor::default_param_value("n_cells_init", 10, 1000), 200);
}

TEST(Options, FillInputOptions) {
    baysor::MoleculeInputOptions molecule_opts;
    molecule_opts.min_molecules_per_cell = 20;
    baysor::fill_and_check_molecule_input_options(molecule_opts);

    baysor::PriorInputOptions prior_opts;
    prior_opts.type = baysor::PriorInputType::Column;
    prior_opts.column_name = "cell_id";
    baysor::fill_and_check_prior_input_options(prior_opts, molecule_opts.min_molecules_per_cell);

    EXPECT_GT(molecule_opts.confidence_nn_id, 0);
    EXPECT_GT(prior_opts.min_molecules_per_segment, 0);
}

TEST(Options, LoadConfigFromToml) {
    auto path = write_temp_csv(
        "[molecules]\n"
        "x = \"pos_x\"\n"
        "min_molecules_per_cell = 30\n"
        "\n"
        "[prior]\n"
        "column_name = \"cell_id\"\n"
        "unassigned_label = \"UNASSIGNED\"\n"
        "\n"
        "[segmentation]\n"
        "scale = 7.5\n"
        "cluster_method = \"louvain\"\n"
        "n_clusters = 6\n",
        ".toml"
    );

    auto opts = baysor::load_config(path);

    EXPECT_EQ(opts.molecules.x_col, "pos_x");
    EXPECT_EQ(opts.molecules.min_molecules_per_cell, 30);
    EXPECT_EQ(opts.prior.column_name, "cell_id");
    EXPECT_EQ(opts.prior.unassigned_label, "UNASSIGNED");
    EXPECT_DOUBLE_EQ(opts.segmentation.scale, 7.5);
    EXPECT_EQ(opts.segmentation.cluster_method, baysor::ClusterMethod::Louvain);
    EXPECT_EQ(opts.segmentation.n_clusters, 6);

    std::remove(path.c_str());
}

TEST(Options, LoadConfigAssignsMethodSpecificClusterDefault) {
    auto path = write_temp_csv(
        "[molecules]\n"
        "min_molecules_per_cell = 30\n"
        "\n"
        "[segmentation]\n"
        "cluster_method = \"louvain\"\n",
        ".toml"
    );

    auto opts = baysor::load_config(path);

    EXPECT_EQ(opts.segmentation.cluster_method, baysor::ClusterMethod::Louvain);
    EXPECT_EQ(opts.segmentation.n_clusters, 10);

    std::remove(path.c_str());
}

TEST(Options, ParseClusterMethodAcceptsLegacyAlias) {
    EXPECT_EQ(baysor::parse_cluster_method("mrf"), baysor::ClusterMethod::Mrf);
    EXPECT_EQ(baysor::parse_cluster_method("ica_mrf"), baysor::ClusterMethod::Mrf);
    EXPECT_EQ(baysor::parse_cluster_method("ica"), baysor::ClusterMethod::Mrf);
    EXPECT_EQ(baysor::parse_cluster_method("leiden"), baysor::ClusterMethod::Leiden);
    EXPECT_EQ(baysor::cluster_method_to_string(baysor::ClusterMethod::Mrf), "mrf");
    EXPECT_EQ(baysor::cluster_method_to_string(baysor::ClusterMethod::Leiden), "leiden");
    EXPECT_EQ(baysor::default_cluster_count(baysor::ClusterMethod::Mrf), 4);
    EXPECT_EQ(baysor::default_cluster_count(baysor::ClusterMethod::Louvain), 10);
    EXPECT_EQ(baysor::default_cluster_count(baysor::ClusterMethod::Leiden), 10);
}

// ============================================================================
// Utility functions
// ============================================================================

TEST(Utils, CountArray) {
    std::vector<int> values = {1, 2, 3, 1, 2, 1};
    auto counts = baysor::count_array(values);

    ASSERT_EQ(counts.size(), 3u);
    EXPECT_EQ(counts[0], 3);  // value 1
    EXPECT_EQ(counts[1], 2);  // value 2
    EXPECT_EQ(counts[2], 1);  // value 3
}

TEST(Utils, SplitIds) {
    std::vector<int> factor = {2, 1, 3, 1, 2};
    auto groups = baysor::split_ids(factor);

    ASSERT_EQ(groups.size(), 3u);
    // Group 1: indices where factor == 1
    ASSERT_EQ(groups[0].size(), 2u);
    EXPECT_EQ(groups[0][0], 1);
    EXPECT_EQ(groups[0][1], 3);
    // Group 2: indices where factor == 2
    ASSERT_EQ(groups[1].size(), 2u);
    EXPECT_EQ(groups[1][0], 0);
    EXPECT_EQ(groups[1][1], 4);
}

TEST(Utils, SplitStringList) {
    auto result = baysor::split_string_list("Blank*, MALAT1 , NegCtrl");
    ASSERT_EQ(result.size(), 3u);
    EXPECT_EQ(result[0], "Blank*");
    EXPECT_EQ(result[1], "MALAT1");
    EXPECT_EQ(result[2], "NegCtrl");
}

TEST(Utils, Strip) {
    EXPECT_EQ(baysor::strip("  hello  "), "hello");
    EXPECT_EQ(baysor::strip("\thello\n"), "hello");
    EXPECT_EQ(baysor::strip(""), "");
}

// ============================================================================
// KNN (nanoflann-based)
// ============================================================================

TEST(KNN, SmallGrid) {
    // 4 points on a unit square
    Eigen::MatrixXd pts(2, 4);
    pts.col(0) << 0.0, 0.0;
    pts.col(1) << 1.0, 0.0;
    pts.col(2) << 0.0, 1.0;
    pts.col(3) << 1.0, 1.0;

    auto result = baysor::knn_parallel(pts, pts, 3, true);

    ASSERT_EQ(result.n, 4);
    ASSERT_EQ(result.k, 3);
    ASSERT_EQ(result.indices.size(), 12u);
    ASSERT_EQ(result.distances.size(), 12u);

    // For point (0,0), nearest is itself (dist=0), then (1,0) or (0,1) (dist=1)
    EXPECT_EQ(result.idx_row(0)[0], 0);  // self
    EXPECT_NEAR(result.dist_row(0)[0], 0.0, 1e-10);
    EXPECT_NEAR(result.dist_row(0)[1], 1.0, 1e-10);  // adjacent
    EXPECT_NEAR(result.dist_row(0)[2], 1.0, 1e-10);  // adjacent
}

TEST(KNN, SelfQuery) {
    Eigen::MatrixXd pts(2, 3);
    pts.col(0) << 0.0, 0.0;
    pts.col(1) << 3.0, 0.0;
    pts.col(2) << 0.0, 4.0;

    auto result = baysor::knn_parallel(pts, pts, 2, true);

    // First neighbor is always self
    for (int i = 0; i < 3; ++i) {
        EXPECT_EQ(result.idx_row(i)[0], i);
        EXPECT_NEAR(result.dist_row(i)[0], 0.0, 1e-10);
    }
}

TEST(KNN, 3D) {
    Eigen::MatrixXd pts(3, 2);
    pts.col(0) << 0.0, 0.0, 0.0;
    pts.col(1) << 3.0, 4.0, 0.0;

    auto result = baysor::knn_parallel(pts, pts, 2, true);

    EXPECT_NEAR(result.dist_row(0)[1], 5.0, 1e-10);
    EXPECT_NEAR(result.dist_row(1)[1], 5.0, 1e-10);
}

TEST(KNN, 3DTieOrderIsDeterministic) {
    Eigen::MatrixXd pts(3, 4);
    pts.col(0) << 0.0, 0.0, 0.0;
    pts.col(1) << 1.0, 0.0, 0.0;
    pts.col(2) << -1.0, 0.0, 0.0;
    pts.col(3) << 0.0, 1.0, 0.0;

    Eigen::MatrixXd query(3, 1);
    query.col(0) = pts.col(0);

    auto result = baysor::knn_parallel(pts, query, 4, true);

    ASSERT_EQ(result.n, 1);
    ASSERT_EQ(result.k, 4);
    EXPECT_EQ(result.idx_row(0)[0], 0);
    EXPECT_EQ(result.idx_row(0)[1], 1);
    EXPECT_EQ(result.idx_row(0)[2], 2);
    EXPECT_EQ(result.idx_row(0)[3], 3);
    EXPECT_NEAR(result.dist_row(0)[1], 1.0, 1e-10);
    EXPECT_NEAR(result.dist_row(0)[2], 1.0, 1e-10);
    EXPECT_NEAR(result.dist_row(0)[3], 1.0, 1e-10);
}

// ============================================================================
// count_array_sparse
// ============================================================================

TEST(SparseCount, Basic) {
    // Gene IDs (1-based): {1, 2, 1, 3, 2}
    int values[] = {1, 2, 1, 3, 2};
    auto sv = baysor::count_array_sparse(values, 5, 3);

    EXPECT_EQ(sv.size(), 3);
    EXPECT_FLOAT_EQ(sv.coeff(0), 2.0f);  // gene 1: 2 counts
    EXPECT_FLOAT_EQ(sv.coeff(1), 2.0f);  // gene 2: 2 counts
    EXPECT_FLOAT_EQ(sv.coeff(2), 1.0f);  // gene 3: 1 count
}

TEST(SparseCount, Weighted) {
    int values[] = {1, 1, 2};
    double weights[] = {0.5, 1.5, 3.0};
    auto sv = baysor::count_array_sparse(values, 3, 2, weights);

    EXPECT_FLOAT_EQ(sv.coeff(0), 2.0f);  // gene 1: 0.5 + 1.5
    EXPECT_FLOAT_EQ(sv.coeff(1), 3.0f);  // gene 2: 3.0
}

TEST(SparseCount, Normalized) {
    int values[] = {1, 2, 2};
    auto sv = baysor::count_array_sparse(values, 3, 2, nullptr, true);

    // gene 1: 1/3, gene 2: 2/3
    EXPECT_NEAR(sv.coeff(0), 1.0f / 3.0f, 1e-5f);
    EXPECT_NEAR(sv.coeff(1), 2.0f / 3.0f, 1e-5f);
}

// ============================================================================
// Triangulation
// ============================================================================

TEST(Triangulation, NormalizePoints) {
    Eigen::MatrixXd pts(2, 4);
    pts.col(0) << 0.0, 0.0;
    pts.col(1) << 100.0, 0.0;
    pts.col(2) << 0.0, 100.0;
    pts.col(3) << 100.0, 100.0;

    auto norm = baysor::normalize_points(pts);

    // All values should be in [1.0, 2.0]
    EXPECT_GE(norm.minCoeff(), 1.0);
    EXPECT_LE(norm.maxCoeff(), 2.0);
}

TEST(Triangulation, AdjacencyList2D) {
    // Create a small 2D grid: 4 points on a unit square
    Eigen::MatrixXd pts(2, 4);
    pts.col(0) << 0.0, 0.0;
    pts.col(1) << 1.0, 0.0;
    pts.col(2) << 0.0, 1.0;
    pts.col(3) << 1.0, 1.0;

    auto result = baysor::adjacency_list(pts, /*filter=*/false);

    // Delaunay of 4 points on a square should produce 5 edges
    // (4 sides + 1 diagonal)
    EXPECT_GE(static_cast<int>(result.edge_src.size()), 4);
    EXPECT_LE(static_cast<int>(result.edge_src.size()), 6);

    // All edges should have positive distance
    for (double d : result.edge_dists) {
        EXPECT_GT(d, 0.0);
    }
}

TEST(Triangulation, AdjacencyList3DKNN) {
    // 3D should use KNN
    Eigen::MatrixXd pts(3, 5);
    pts.col(0) << 0.0, 0.0, 0.0;
    pts.col(1) << 1.0, 0.0, 0.0;
    pts.col(2) << 0.0, 1.0, 0.0;
    pts.col(3) << 0.0, 0.0, 1.0;
    pts.col(4) << 1.0, 1.0, 1.0;

    auto result = baysor::adjacency_list(pts, /*filter=*/false);

    EXPECT_GT(static_cast<int>(result.edge_src.size()), 0);
}

// ============================================================================
// Boundary estimation
// ============================================================================

static std::vector<std::pair<double, double>> canonical_polygon_vertices(const Eigen::MatrixXd& poly) {
    std::vector<std::pair<double, double>> vertices;
    vertices.reserve(static_cast<size_t>(poly.cols()));
    for (int i = 0; i < poly.cols(); ++i) {
        vertices.push_back({poly(0, i), poly(1, i)});
    }
    if (vertices.empty()) return vertices;

    auto canonicalize_direction = [](std::vector<std::pair<double, double>> vals) {
        const auto it = std::min_element(vals.begin(), vals.end());
        std::rotate(vals.begin(), it, vals.end());
        return vals;
    };

    auto forward = canonicalize_direction(vertices);
    std::reverse(vertices.begin(), vertices.end());
    auto backward = canonicalize_direction(vertices);
    return std::lexicographical_compare(backward.begin(), backward.end(), forward.begin(), forward.end())
        ? backward
        : forward;
}

static void expect_polygon_collection_near(
    const baysor::PolygonCollection& actual,
    const baysor::PolygonCollection& expected,
    double tol = 1e-10
) {
    ASSERT_EQ(actual.size(), expected.size());
    for (const auto& [cell, expected_poly] : expected) {
        auto it = actual.find(cell);
        ASSERT_NE(it, actual.end()) << cell;
        auto actual_vertices = canonical_polygon_vertices(it->second);
        auto expected_vertices = canonical_polygon_vertices(expected_poly);
        ASSERT_EQ(actual_vertices.size(), expected_vertices.size()) << cell;
        for (size_t i = 0; i < expected_vertices.size(); ++i) {
            EXPECT_NEAR(actual_vertices[i].first, expected_vertices[i].first, tol) << cell;
            EXPECT_NEAR(actual_vertices[i].second, expected_vertices[i].second, tol) << cell;
        }
    }
}

TEST(BoundaryEstimation, BoundaryPolygonsAuto3DPerZ) {
    Eigen::MatrixXd pos(3, 8);
    pos.col(0) << 0.0, 0.0, 0.0;
    pos.col(1) << 2.0, 0.0, 0.0;
    pos.col(2) << 2.0, 2.0, 0.0;
    pos.col(3) << 0.0, 2.0, 0.0;
    pos.col(4) << 10.0, 10.0, 1.0;
    pos.col(5) << 12.0, 10.0, 1.0;
    pos.col(6) << 12.0, 12.0, 1.0;
    pos.col(7) << 10.0, 12.0, 1.0;

    std::vector<int> assignment = {1, 1, 1, 1, 2, 2, 2, 2};
    std::vector<std::string> cell_names = {"cell_1", "cell_2"};

    auto [joined, stack] = baysor::boundary_polygons_auto(
        pos, assignment, /*estimate_per_z=*/true, &cell_names, /*verbose=*/false);

    ASSERT_EQ(joined.size(), 2u);
    EXPECT_GT(joined.at("cell_1").cols(), 0);
    EXPECT_GT(joined.at("cell_2").cols(), 0);

    ASSERT_GE(stack.size(), 3u);
    EXPECT_EQ(stack[0].first, "2d");

    bool saw_z0 = false;
    bool saw_z1 = false;
    for (const auto& [layer, polys] : stack) {
        if (layer == "0") {
            saw_z0 = true;
            ASSERT_EQ(polys.size(), 1u);
            EXPECT_TRUE(polys.count("cell_1"));
        } else if (layer == "1") {
            saw_z1 = true;
            ASSERT_EQ(polys.size(), 1u);
            EXPECT_TRUE(polys.count("cell_2"));
        }
    }
    EXPECT_TRUE(saw_z0);
    EXPECT_TRUE(saw_z1);
}

TEST(BoundaryEstimation, BoundaryPolygonsAutoStableAcrossThreadCounts) {
    constexpr int n_cells = 4;
    constexpr int points_per_cell = 4;
    Eigen::MatrixXd pos(3, n_cells * points_per_cell);
    std::vector<int> assignment(n_cells * points_per_cell);
    std::vector<std::string> cell_names = {"cell_1", "cell_2", "cell_3", "cell_4"};

    auto add_square = [&](int cell, double x0, double y0, double z) {
        const int offset = (cell - 1) * points_per_cell;
        pos.col(offset + 0) << x0, y0, z;
        pos.col(offset + 1) << x0 + 2.0, y0, z;
        pos.col(offset + 2) << x0 + 2.0, y0 + 2.0, z;
        pos.col(offset + 3) << x0, y0 + 2.0, z;
        for (int i = 0; i < points_per_cell; ++i) assignment[offset + i] = cell;
    };

    add_square(1, 0.0, 0.0, 0.0);
    add_square(2, 10.0, 0.0, 0.0);
    add_square(3, 0.0, 10.0, 1.0);
    add_square(4, 10.0, 10.0, 1.0);

    int old_threads = baysor::thread_pool_size();
    baysor::set_thread_pool_size(1);
    auto [joined_1, stack_1] = baysor::boundary_polygons_auto(
        pos, assignment, /*estimate_per_z=*/true, &cell_names, /*verbose=*/false);
    baysor::set_thread_pool_size(4);
    auto [joined_4, stack_4] = baysor::boundary_polygons_auto(
        pos, assignment, /*estimate_per_z=*/true, &cell_names, /*verbose=*/false);
    baysor::set_thread_pool_size(old_threads);

    expect_polygon_collection_near(joined_4, joined_1);
    ASSERT_EQ(stack_4.size(), stack_1.size());
    for (size_t i = 0; i < stack_1.size(); ++i) {
        EXPECT_EQ(stack_4[i].first, stack_1[i].first);
        expect_polygon_collection_near(stack_4[i].second, stack_1[i].second);
    }
}

TEST(BoundaryEstimation, BoundaryPolygonsFromGridKeepsTouchingLabelsSeparate) {
    Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> grid(3, 4);
    grid <<
        1, 1, 2, 2,
        1, 1, 2, 2,
        0, 0, 0, 0;

    auto polys = baysor::boundary_polygons_from_grid(grid);

    ASSERT_GE(polys.size(), 2u);
    EXPECT_GT(polys[0].cols(), 0);
    EXPECT_GT(polys[1].cols(), 0);
}

TEST(Output, SavePolygonsGeoJsonKeepsTriangleFeature) {
    baysor::PolygonCollection polygons;
    Eigen::MatrixXd tri(2, 3);
    tri <<
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
    polygons["cell_1"] = tri;

    const std::string path = write_temp_csv("", ".json");
    baysor::save_polygons_geojson(polygons, path, "FeatureCollection");

    std::ifstream f(path);
    ASSERT_TRUE(f.good());
    std::stringstream buf;
    buf << f.rdbuf();
    const std::string s = buf.str();

    EXPECT_NE(s.find("\"type\":\"FeatureCollection\""), std::string::npos);
    EXPECT_NE(s.find("\"id\":\"cell_1\""), std::string::npos);
    EXPECT_NE(s.find("\"properties\":{\"cell\":\"cell_1\"}"), std::string::npos);
}

TEST(Output, ParseOutputStyle) {
    EXPECT_EQ(baysor::parse_output_style("legacy"), baysor::OutputStyle::Legacy);
    EXPECT_EQ(baysor::parse_output_style("parquet"), baysor::OutputStyle::Parquet);
    EXPECT_THROW(baysor::parse_output_style("xenium"), std::invalid_argument);
    EXPECT_THROW(baysor::parse_output_style("csv"), std::invalid_argument);
}

TEST(Output, GetOutputPathsSupportsLegacyAndParquetStyles) {
    auto legacy = baysor::get_output_paths("out_dir/", baysor::OutputStyle::Legacy, "tsv");
    EXPECT_EQ(legacy.segmented_df, "out_dir/segmentation.csv");
    EXPECT_EQ(legacy.cell_stats, "out_dir/segmentation_cell_stats.csv");
    EXPECT_EQ(legacy.polygons_2d, "out_dir/segmentation_polygons_2d.json");
    EXPECT_EQ(legacy.polygons_3d, "out_dir/segmentation_polygons_3d.json");
    EXPECT_EQ(legacy.counts, "out_dir/segmentation_counts.tsv");
    EXPECT_EQ(legacy.params_dump, "out_dir/segmentation_params.dump.toml");
    EXPECT_EQ(legacy.log_file, "out_dir/segmentation_log.log");

    auto parquet = baysor::get_output_paths("out_dir/", baysor::OutputStyle::Parquet, "loom");
    EXPECT_EQ(parquet.segmented_df, "out_dir/molecules.parquet");
    EXPECT_EQ(parquet.cell_stats, "out_dir/cells.parquet");
    EXPECT_EQ(parquet.polygons_2d, "out_dir/cell_boundaries.parquet");
    EXPECT_EQ(parquet.polygons_3d, "out_dir/cell_boundaries_3d.parquet");
    EXPECT_EQ(parquet.counts, "out_dir/feature_matrix.h5");
    EXPECT_EQ(parquet.params_dump, "out_dir/run_params.toml");
    EXPECT_EQ(parquet.log_file, "out_dir/run.log");

}

TEST(Output, SaveSegmentedDfAddsTranscriptIdForXeniumLikeInput) {
    baysor::MoleculeData data;
    data.x = {1.0, 2.0};
    data.y = {3.0, 4.0};
    data.gene = {1, 1};
    data.gene_names = {"MALAT1"};
    data.source_transcript_id = {281474976710656ULL, 281474976710657ULL};

    std::vector<int> assignment = {1, 0};
    const std::string path = write_temp_csv("", ".csv");
    baysor::save_segmented_df(data, assignment, data.gene_names, path);

    std::ifstream f(path);
    ASSERT_TRUE(f.good());
    std::stringstream buf;
    buf << f.rdbuf();
    const std::string s = buf.str();

    EXPECT_NE(s.find("transcript_id,cell,gene,x,y,is_noise"), std::string::npos);
    EXPECT_NE(s.find("281474976710656,cell_1,MALAT1,1,3,false"), std::string::npos);
    EXPECT_NE(s.find("281474976710657,0,MALAT1,2,4,true"), std::string::npos);
}

TEST(Output, SavePolygonsGeoParquetWritesGeometryMetadata) {
    baysor::PolygonCollection polygons;
    Eigen::MatrixXd tri(2, 3);
    tri <<
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
    polygons["cell_1"] = tri;

    const std::string path = write_temp_csv("", ".parquet");
    baysor::save_polygons_geoparquet(polygons, path);

    auto infile = arrow_test_unwrap(arrow::io::ReadableFile::Open(path));
    auto reader = arrow_test_unwrap(parquet::arrow::OpenFile(infile, arrow::default_memory_pool()));
    auto pq_reader = parquet::ParquetFileReader::Open(infile);
    auto kv = pq_reader->metadata()->key_value_metadata();
    ASSERT_NE(kv, nullptr);
    ASSERT_GE(kv->size(), 1);
    bool saw_geo = false;
    for (int i = 0; i < kv->size(); ++i) {
        if (kv->key(i) == "geo") {
            saw_geo = true;
            EXPECT_NE(kv->value(i).find("\"primary_column\":\"geometry\""), std::string::npos);
        }
    }
    EXPECT_TRUE(saw_geo);

    std::shared_ptr<arrow::Table> table;
    arrow_expect_ok(reader->ReadTable(&table));
    ASSERT_NE(table, nullptr);
    EXPECT_EQ(table->num_rows(), 1);
    EXPECT_EQ(table->num_columns(), 3);
    EXPECT_EQ(table->schema()->field(0)->name(), "cell");
    EXPECT_EQ(table->schema()->field(1)->name(), "n_vertices");
    EXPECT_EQ(table->schema()->field(2)->name(), "geometry");
}

TEST(Output, SaveMatrixTo10xH5WritesExpectedStructure) {
    Eigen::SparseMatrix<float> matrix(2, 3);  // cells x genes
    std::vector<Eigen::Triplet<float>> trips = {
        {0, 0, 2.0f},
        {0, 2, 5.0f},
        {1, 1, 7.0f}
    };
    matrix.setFromTriplets(trips.begin(), trips.end());
    matrix.makeCompressed();

    const std::vector<std::string> gene_names = {"G1", "G2", "G3"};
    const std::vector<std::string> cell_names = {"cell_1", "cell_2"};
    const std::string path = write_temp_csv("", ".h5");

    baysor::save_matrix_to_10x_h5(matrix, gene_names, cell_names, path);

    hid_t fid = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(fid, 0);

    hid_t shape_ds = H5Dopen2(fid, "/matrix/shape", H5P_DEFAULT);
    ASSERT_GE(shape_ds, 0);
    std::vector<int64_t> shape(2, -1);
    ASSERT_GE(H5Dread(shape_ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, shape.data()), 0);
    EXPECT_EQ(shape[0], 3);
    EXPECT_EQ(shape[1], 2);
    H5Dclose(shape_ds);

    hid_t data_ds = H5Dopen2(fid, "/matrix/data", H5P_DEFAULT);
    ASSERT_GE(data_ds, 0);
    hid_t data_space = H5Dget_space(data_ds);
    hsize_t data_dims[1] = {0};
    ASSERT_EQ(H5Sget_simple_extent_ndims(data_space), 1);
    ASSERT_EQ(H5Sget_simple_extent_dims(data_space, data_dims, nullptr), 1);
    EXPECT_EQ(data_dims[0], 3);
    H5Sclose(data_space);
    H5Dclose(data_ds);

    hid_t feat_name_ds = H5Dopen2(fid, "/matrix/features/name", H5P_DEFAULT);
    ASSERT_GE(feat_name_ds, 0);
    hid_t feat_name_space = H5Dget_space(feat_name_ds);
    hsize_t feat_dims[1] = {0};
    ASSERT_EQ(H5Sget_simple_extent_ndims(feat_name_space), 1);
    ASSERT_EQ(H5Sget_simple_extent_dims(feat_name_space, feat_dims, nullptr), 1);
    EXPECT_EQ(feat_dims[0], 3);
    H5Sclose(feat_name_space);
    H5Dclose(feat_name_ds);

    H5Fclose(fid);
}

TEST(Output, SaveMatrixToLoomWritesGeneByCellMatrix) {
    Eigen::SparseMatrix<float, Eigen::RowMajor> matrix(2, 3);  // cells x genes
    std::vector<Eigen::Triplet<float>> trips = {
        {0, 0, 2.0f},
        {0, 2, 5.0f},
        {1, 1, 7.0f}
    };
    matrix.setFromTriplets(trips.begin(), trips.end());
    matrix.makeCompressed();

    const std::vector<std::string> gene_names = {"G1", "G2", "G3"};
    const std::vector<std::string> cell_names = {"cell_1", "cell_2"};
    const std::string path = write_temp_csv("", ".loom");

    baysor::save_matrix_to_loom(matrix, gene_names, cell_names, path);

    hid_t fid = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(fid, 0);

    hid_t matrix_ds = H5Dopen2(fid, "/matrix", H5P_DEFAULT);
    ASSERT_GE(matrix_ds, 0);
    hid_t matrix_space = H5Dget_space(matrix_ds);
    hsize_t matrix_dims[2] = {0, 0};
    ASSERT_EQ(H5Sget_simple_extent_ndims(matrix_space), 2);
    ASSERT_EQ(H5Sget_simple_extent_dims(matrix_space, matrix_dims, nullptr), 2);
    EXPECT_EQ(matrix_dims[0], 3);
    EXPECT_EQ(matrix_dims[1], 2);

    std::vector<float> values(6, -1.0f);
    ASSERT_GE(H5Dread(matrix_ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()), 0);
    const std::vector<float> expected = {
        2.0f, 0.0f,
        0.0f, 7.0f,
        5.0f, 0.0f
    };
    EXPECT_EQ(values, expected);
    H5Sclose(matrix_space);
    H5Dclose(matrix_ds);

    hid_t row_name_ds = H5Dopen2(fid, "/row_attrs/Name", H5P_DEFAULT);
    ASSERT_GE(row_name_ds, 0);
    hid_t row_name_space = H5Dget_space(row_name_ds);
    hsize_t row_dims[1] = {0};
    ASSERT_EQ(H5Sget_simple_extent_ndims(row_name_space), 1);
    ASSERT_EQ(H5Sget_simple_extent_dims(row_name_space, row_dims, nullptr), 1);
    EXPECT_EQ(row_dims[0], 3);
    H5Sclose(row_name_space);
    H5Dclose(row_name_ds);

    hid_t col_name_ds = H5Dopen2(fid, "/col_attrs/Name", H5P_DEFAULT);
    ASSERT_GE(col_name_ds, 0);
    hid_t col_name_space = H5Dget_space(col_name_ds);
    hsize_t col_dims[1] = {0};
    ASSERT_EQ(H5Sget_simple_extent_ndims(col_name_space), 1);
    ASSERT_EQ(H5Sget_simple_extent_dims(col_name_space, col_dims, nullptr), 1);
    EXPECT_EQ(col_dims[0], 2);
    H5Sclose(col_name_space);
    H5Dclose(col_name_ds);

    hid_t cell_id_ds = H5Dopen2(fid, "/col_attrs/CellID", H5P_DEFAULT);
    ASSERT_GE(cell_id_ds, 0);
    hid_t cell_id_space = H5Dget_space(cell_id_ds);
    hsize_t cell_id_dims[1] = {0};
    ASSERT_EQ(H5Sget_simple_extent_ndims(cell_id_space), 1);
    ASSERT_EQ(H5Sget_simple_extent_dims(cell_id_space, cell_id_dims, nullptr), 1);
    EXPECT_EQ(cell_id_dims[0], 2);
    H5Sclose(cell_id_space);
    H5Dclose(cell_id_ds);

    H5Fclose(fid);
}

TEST(Output, SaveMatrixToLoomRejectsMismatchedColAttrLength) {
    Eigen::SparseMatrix<float> matrix(2, 3);  // cells x genes
    const std::vector<std::string> gene_names = {"G1", "G2", "G3"};
    const std::vector<std::string> cell_names = {"cell_1", "cell_2"};
    const std::string path = write_temp_csv("", ".loom");

    baysor::LoomColAttrs col_attrs;
    col_attrs["confidence"] = std::vector<double>{0.5};

    EXPECT_THROW(
        baysor::save_matrix_to_loom(matrix, gene_names, cell_names, path, col_attrs),
        std::runtime_error);
}

TEST(Xenium, ManifestHelpers) {
    EXPECT_TRUE(baysor::is_xenium_manifest_path("experiment.xenium"));
    EXPECT_FALSE(baysor::is_xenium_manifest_path("transcripts.parquet"));

    char tmpl[] = "/tmp/baysor_xenium_manifest_XXXXXX";
    char* created = mkdtemp(tmpl);
    ASSERT_NE(created, nullptr);

    const std::filesystem::path dataset_dir(created);
    const std::string manifest_path = (dataset_dir / "experiment.xenium").string();
    {
        std::ofstream manifest(manifest_path);
        ASSERT_TRUE(static_cast<bool>(manifest));
        manifest << "{}\n";
    }
    {
        std::ofstream transcripts((dataset_dir / "transcripts.parquet").string(), std::ios::binary);
        ASSERT_TRUE(static_cast<bool>(transcripts));
    }

    auto ctx = baysor::load_xenium_manifest_context(manifest_path);
    EXPECT_EQ(ctx.manifest_path, manifest_path);
    EXPECT_EQ(ctx.dataset_dir, dataset_dir.string());
    EXPECT_EQ(ctx.transcripts_path, (dataset_dir / "transcripts.parquet").string());

    std::error_code ec;
    std::filesystem::remove_all(dataset_dir, ec);
    EXPECT_FALSE(ec) << ec.message();
}

// ============================================================================
// AdjList::from_edge_list
// ============================================================================

TEST(AdjList, FromEdgeList) {
    // 3 vertices, 2 edges: 0-1, 1-2
    int src[] = {0, 1};
    int dst[] = {1, 2};
    double wts[] = {0.5, 0.8};

    auto adj = baysor::AdjList::from_edge_list(src, dst, wts, 2, 3);

    EXPECT_EQ(adj.n_molecules(), 3);

    // Vertex 0: neighbors should include 1
    EXPECT_GE(adj.neighbor_count(0), 1);
    // Vertex 1: neighbors should include 0 and 2
    EXPECT_GE(adj.neighbor_count(1), 2);
    // Vertex 2: neighbors should include 1
    EXPECT_GE(adj.neighbor_count(2), 1);
}

TEST(AdjList, EmptyGraph) {
    auto adj = baysor::AdjList::from_edge_list(nullptr, nullptr, nullptr, 0, 5);

    EXPECT_EQ(adj.n_molecules(), 5);
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(adj.neighbor_count(i), 0);
    }
}

// ============================================================================
// build_molecule_graph
// ============================================================================

TEST(MoleculeGraph, Basic) {
    // Use a 2D grid to avoid degenerate collinear triangulation
    baysor::MoleculeData data;
    // 3x3 grid
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            data.x.push_back(i * 1.0);
            data.y.push_back(j * 1.0);
            data.gene.push_back(1);
        }
    }

    auto adj = baysor::build_molecule_graph(data);

    EXPECT_EQ(adj.n_molecules(), 9);
    // Center point (index 4) should have neighbors
    EXPECT_GE(adj.neighbor_count(4), 2);
    // Total edges should be reasonable
    EXPECT_GT(adj.nnz(), 0);
}

// ============================================================================
// Noise estimation
// ============================================================================

TEST(NoiseEstimation, FitBimodal) {
    // Create synthetic bimodal data: signal (low distances) + noise (high distances)
    std::vector<double> edge_lengths;
    std::mt19937 rng(42);

    // Signal: Normal(2.0, 0.5) — 100 points
    std::normal_distribution<double> signal_dist(2.0, 0.5);
    for (int i = 0; i < 100; ++i) {
        edge_lengths.push_back(std::max(0.1, signal_dist(rng)));
    }
    // Noise: Normal(8.0, 1.0) — 50 points
    std::normal_distribution<double> noise_dist(8.0, 1.0);
    for (int i = 0; i < 50; ++i) {
        edge_lengths.push_back(noise_dist(rng));
    }

    int n = static_cast<int>(edge_lengths.size());

    // Build a simple fully-connected adjacency for testing
    // (each point connected to all others — won't happen in practice but tests the EM)
    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n && j < i + 5; ++j) {
            src.push_back(i);
            dst.push_back(j);
            wts.push_back(1.0);
        }
    }
    auto adj = baysor::AdjList::from_edge_list(src.data(), dst.data(), wts.data(),
                                                static_cast<int>(src.size()), n);

    auto result = baysor::fit_noise_probabilities(edge_lengths, adj, nullptr, 100, 0.005, false);

    // Signal mu should be close to 2.0, noise mu close to 8.0
    EXPECT_NEAR(result.signal_mu, 2.0, 1.5);
    EXPECT_NEAR(result.noise_mu, 8.0, 2.0);

    // Assignment probs should be n x 2
    EXPECT_EQ(result.assignment_probs.rows(), n);
    EXPECT_EQ(result.assignment_probs.cols(), 2);

    // Signal points should have high signal probability
    double avg_signal_prob = 0.0;
    for (int i = 0; i < 100; ++i) {
        avg_signal_prob += result.assignment_probs(i, 0);
    }
    avg_signal_prob /= 100.0;
    EXPECT_GT(avg_signal_prob, 0.5);

    // Noise points should have low signal probability
    double avg_noise_prob = 0.0;
    for (int i = 100; i < n; ++i) {
        avg_noise_prob += result.assignment_probs(i, 0);
    }
    avg_noise_prob /= 50.0;
    EXPECT_LT(avg_noise_prob, 0.5);
}

TEST(NoiseEstimation, AppendConfidence) {
    // Create a grid of molecules with varied spacing
    baysor::MoleculeData data;
    // Dense cluster (signal-like)
    for (int i = 0; i < 20; ++i) {
        data.x.push_back(i * 0.5);
        data.y.push_back(0.0);
        data.gene.push_back(1);
    }
    // Sparse points (noise-like)
    for (int i = 0; i < 10; ++i) {
        data.x.push_back(100.0 + i * 10.0);
        data.y.push_back(100.0);
        data.gene.push_back(1);
    }
    data.gene_names = {"GeneA"};

    baysor::append_confidence(data, 3);

    ASSERT_EQ(data.confidence.size(), 30u);

    // All confidence values should be in [0, 1]
    for (double c : data.confidence) {
        EXPECT_GE(c, 0.0);
        EXPECT_LE(c, 1.0);
    }

    // Dense cluster should generally have higher confidence than sparse points
    double avg_dense = 0.0;
    for (int i = 0; i < 20; ++i) avg_dense += data.confidence[i];
    avg_dense /= 20.0;

    double avg_sparse = 0.0;
    for (int i = 20; i < 30; ++i) avg_sparse += data.confidence[i];
    avg_sparse /= 10.0;

    EXPECT_GT(avg_dense, avg_sparse);
}

TEST(NoiseEstimation, JuliaParityWithoutPrior) {
    const std::vector<double> edge_lengths = {
        0.08499664857992582, 0.7137346903406472, 0.3061761980837148,
        0.7206726148556293,  0.09408790764541042, 0.06353973485938673,
        0.12871716807062483, 0.3445578901832464,  0.3269150839207085,
        0.6461584760429055
    };
    const std::vector<double> expected_signal_probs = {
        1.0, 4.564161337603029e-5, 1.0, 4.060842872158915e-5, 1.0,
        1.0, 1.0, 1.0, 1.0, 0.0011055181085140892
    };
    const std::vector<int> expected_assignment = {1, 2, 1, 2, 1, 1, 1, 1, 1, 2};

    auto adj = make_chain_adj_list(static_cast<int>(edge_lengths.size()));
    auto result = baysor::fit_noise_probabilities(
        edge_lengths, adj, nullptr, /*max_iters=*/10000, /*tol=*/0.005, /*verbose=*/false);

    std::vector<double> signal_probs(edge_lengths.size());
    for (int i = 0; i < result.assignment_probs.rows(); ++i) {
        signal_probs[i] = result.assignment_probs(i, 0);
    }

    expect_vector_near(signal_probs, expected_signal_probs, 1e-6);
    EXPECT_EQ(result.assignment, expected_assignment);
    EXPECT_NEAR(result.signal_mu, 0.19306303931238833, 1e-6);
    EXPECT_NEAR(result.signal_sigma, 0.1177901622243064, 1e-6);
    EXPECT_NEAR(result.noise_mu, 0.6935851572866298, 1e-6);
    EXPECT_NEAR(result.noise_sigma, 0.03358834820916101, 1e-6);
}

TEST(NoiseEstimation, JuliaParityWithPriorFloor) {
    const std::vector<double> edge_lengths = {
        0.3124500990292511,  0.3050153456358032,  0.2918760634633287,
        0.603571308954251,   0.07861479828819232, 0.08839793695093187,
        0.07366550362329821, 0.05056330797820069, 0.7831762164017797,
        0.591625542529262
    };
    const std::vector<double> min_confidence = {
        0.0, 0.25, 0.0, 0.25, 0.25, 0.04, 0.0, 0.04, 0.0, 0.25
    };
    const std::vector<double> expected_signal_probs = {
        1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 0.0, 1.0
    };
    const std::vector<int> expected_assignment = {1, 1, 1, 1, 1, 1, 1, 1, 2, 1};

    auto adj = make_chain_adj_list(static_cast<int>(edge_lengths.size()));
    auto result = baysor::fit_noise_probabilities(
        edge_lengths, adj, &min_confidence, /*max_iters=*/10000, /*tol=*/0.005, /*verbose=*/false);

    std::vector<double> signal_probs(edge_lengths.size());
    for (int i = 0; i < result.assignment_probs.rows(); ++i) {
        signal_probs[i] = result.assignment_probs(i, 0);
    }

    expect_vector_near(signal_probs, expected_signal_probs, 1e-6);
    EXPECT_EQ(result.assignment, expected_assignment);
    EXPECT_NEAR(result.signal_mu, 0.2661977673836133, 1e-6);
    EXPECT_NEAR(result.signal_sigma, 0.2039598662891038, 1e-6);
    EXPECT_NEAR(result.noise_mu, 0.7831762164017797, 1e-6);
    EXPECT_LE(result.noise_sigma, 1e-9);
}

// ============================================================================
// Molecule clustering
// ============================================================================

TEST(MoleculeClustering, JuliaParityOnMrfWithExplicitInit) {
    const std::vector<int> genes = {1, 1, 1, 2, 2, 3, 3, 3};
    const std::vector<double> confidence = {1.0, 0.9, 0.95, 0.8, 0.85, 0.9, 1.0, 0.95};

    const int edge_src[] = {0, 1, 2, 3, 4, 5, 6};
    const int edge_dst[] = {1, 2, 3, 4, 5, 6, 7};
    const double edge_wt[] = {1.0, 1.2, 0.8, 1.0, 0.8, 1.1, 1.0};
    auto adj = baysor::AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 7, 8);

    Eigen::MatrixXd exprs_init(2, 3);
    exprs_init <<
        0.72, 0.18, 0.10,
        0.08, 0.22, 0.70;

    Eigen::MatrixXd expected_exprs(2, 3);
    expected_exprs <<
        0.6880488047726171, 0.22611659637212975, 0.08583459885525314,
        0.0872518752183331, 0.22286186039622793, 0.689886264385439;

    Eigen::MatrixXd expected_probs(2, 8);
    expected_probs <<
        0.8524769440922572, 0.9286651962470087, 0.8854584870743343, 0.6067947028272389,
        0.4068617974654377, 0.11519259661126363, 0.08034346639566109, 0.13845826934945468,
        0.14752305590774284, 0.0713348037529913, 0.11454151292566576, 0.3932052971727611,
        0.5931382025345622, 0.8848074033887364, 0.9196565336043389, 0.8615417306505453;

    const std::vector<double> expected_diffs = {
        0.1227218878506428, 0.05721821168580932, 0.017336833800294008,
        0.012883262755955301, 0.005312197578650532, 0.005067105452968307,
        0.002731189991373914, 0.002258770150128969, 0.0013986585304166521,
        0.0010752442593115667, 0.000715879234751518, 0.0005307462395459428,
        0.00036679319098976146, 0.00026724685155398716, 0.00018822481673096969,
        0.00013603357451002495, 9.673100527441969e-5, 6.965237532695345e-5,
        4.976818471141842e-5, 3.577802568009991e-5, 2.5626698769198874e-5,
        1.841012471356096e-5, 1.3202986708307774e-5, 9.482323190662666e-6,
        6.804656016851096e-6
    };
    const std::vector<int> expected_assignment = {1, 1, 1, 1, 2, 2, 2, 2};

    auto result = baysor::cluster_molecules_on_mrf(
        genes, adj, confidence, /*n_clusters=*/2,
        /*tol=*/0.0, /*mrf_weight=*/1.0, /*max_iters=*/25, /*verbose=*/false, &exprs_init);

    EXPECT_EQ(result.assignment, expected_assignment);
    expect_matrix_near(result.exprs, expected_exprs, 1e-6);
    expect_matrix_near(result.assignment_probs, expected_probs, 1e-6);
    expect_vector_near(result.diffs, expected_diffs, 1e-6);
    ASSERT_EQ(result.change_fracs.size(), 25u);
    for (double cf : result.change_fracs) {
        EXPECT_DOUBLE_EQ(cf, 1.0);
    }
}

TEST(RandomUtils, XoshiroMatchesJuliaSeed1) {
    baysor::Xoshiro256pp rng(1);
    const std::vector<double> expected = {
        0.07336635446929285,
        0.34924148955718615,
        0.6988266836914685,
        0.6282647403425017,
        0.9149290036628314,
    };

    for (double ex : expected) {
        EXPECT_NEAR(rng.rand_float64(), ex, 1e-15);
    }
}

TEST(MoleculeClustering, IcaWrapperMatchesJuliaPartition) {
    const std::vector<int> genes = {1, 2, 3, 3, 4, 1, 4, 4, 3, 1, 3, 2};
    const std::vector<double> confidence(12, 1.0);

    const int edge_src[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    const int edge_dst[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    const double edge_wt[] = {1.0, 1.3, 1.6, 1.9, 1.0, 1.3, 1.6, 1.9, 1.0, 1.3, 1.6};
    auto adj = baysor::AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 11, 12);

    const std::vector<int> expected_assignment = {2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2};

    auto result = baysor::cluster_molecules_ica(
        genes, adj, confidence, /*n_clusters=*/2,
        /*tol=*/0.0, /*mrf_weight=*/1.0, /*max_iters=*/25, /*verbose=*/false);

    auto swapped = swap_binary_labels(result.assignment);
    EXPECT_TRUE(result.assignment == expected_assignment || swapped == expected_assignment)
        << "actual=" << format_int_vector(result.assignment)
        << " swapped=" << format_int_vector(swapped);
}

TEST(MoleculeClustering, LouvainPartitionSeparatesWeaklyConnectedCliques) {
    const int edge_src[] = {0, 1, 0, 3, 4, 3, 2};
    const int edge_dst[] = {1, 2, 2, 4, 5, 5, 3};
    const double edge_wt[] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.05};
    auto adj = baysor::AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 7, 6);

    auto membership = baysor::louvain_partition(adj, /*resolution=*/1.0, /*max_passes=*/100);

    ASSERT_EQ(membership.size(), 6u);
    EXPECT_EQ(membership[0], membership[1]);
    EXPECT_EQ(membership[1], membership[2]);
    EXPECT_EQ(membership[3], membership[4]);
    EXPECT_EQ(membership[4], membership[5]);
    EXPECT_NE(membership[2], membership[3]);
}

TEST(MoleculeClustering, LouvainPartitionKeepsCompleteGraphTogether) {
    constexpr int n = 12;
    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            src.push_back(i);
            dst.push_back(j);
            wts.push_back(1.0);
        }
    }

    auto adj = baysor::AdjList::from_edge_list(
        src.data(), dst.data(), wts.data(), static_cast<int>(src.size()), n);
    auto membership = baysor::louvain_partition(adj, /*resolution=*/1.0, /*max_passes=*/100);

    ASSERT_EQ(membership.size(), static_cast<size_t>(n));
    for (int i = 1; i < n; ++i) {
        EXPECT_EQ(membership[i], membership[0]);
    }
}

TEST(MoleculeClustering, LouvainPartitionFindsThreeDenseBlocks) {
    constexpr int blocks = 3;
    constexpr int block_size = 8;
    constexpr int n = blocks * block_size;
    std::vector<int> src, dst;
    std::vector<double> wts;

    for (int b = 0; b < blocks; ++b) {
        const int start = b * block_size;
        const int end = start + block_size;
        for (int i = start; i < end; ++i) {
            for (int j = i + 1; j < end; ++j) {
                src.push_back(i);
                dst.push_back(j);
                wts.push_back(1.0);
            }
        }
    }
    // Weak inter-block links.
    for (int b = 0; b < blocks; ++b) {
        const int b2 = (b + 1) % blocks;
        src.push_back(b * block_size);
        dst.push_back(b2 * block_size);
        wts.push_back(0.02);
    }

    auto adj = baysor::AdjList::from_edge_list(
        src.data(), dst.data(), wts.data(), static_cast<int>(src.size()), n);
    auto membership = baysor::louvain_partition(adj, /*resolution=*/1.0, /*max_passes=*/100);

    ASSERT_EQ(membership.size(), static_cast<size_t>(n));
    std::set<int> groups(membership.begin(), membership.end());
    EXPECT_EQ(groups.size(), 3u);
    for (int b = 0; b < blocks; ++b) {
        const int start = b * block_size;
        for (int i = start + 1; i < start + block_size; ++i) {
            EXPECT_EQ(membership[i], membership[start]);
        }
    }
}

TEST(MoleculeClustering, LouvainResolutionSplitsHierarchicalBlocks) {
    constexpr int superblocks = 3;
    constexpr int subblocks_per_super = 2;
    constexpr int block_size = 6;
    constexpr int blocks = superblocks * subblocks_per_super;
    constexpr int n = blocks * block_size;
    std::vector<int> src, dst;
    std::vector<double> wts;

    auto add_clique = [&](int start, int size, double w) {
        for (int i = start; i < start + size; ++i) {
            for (int j = i + 1; j < start + size; ++j) {
                src.push_back(i);
                dst.push_back(j);
                wts.push_back(w);
            }
        }
    };
    auto add_complete_bipartite = [&](int start_a, int start_b, int size, double w) {
        for (int i = start_a; i < start_a + size; ++i) {
            for (int j = start_b; j < start_b + size; ++j) {
                src.push_back(i);
                dst.push_back(j);
                wts.push_back(w);
            }
        }
    };

    for (int b = 0; b < blocks; ++b) {
        add_clique(b * block_size, block_size, 1.0);
    }
    for (int s = 0; s < superblocks; ++s) {
        const int a = (2 * s + 0) * block_size;
        const int b = (2 * s + 1) * block_size;
        add_complete_bipartite(a, b, block_size, 0.35);
    }
    // Weakly connect the three superblocks in a chain.
    for (int s = 0; s < superblocks - 1; ++s) {
        src.push_back((2 * s) * block_size);
        dst.push_back((2 * (s + 1)) * block_size);
        wts.push_back(0.01);
    }

    auto adj = baysor::AdjList::from_edge_list(
        src.data(), dst.data(), wts.data(), static_cast<int>(src.size()), n);
    auto low_res = baysor::louvain_partition(adj, /*resolution=*/0.25, /*max_passes=*/100);
    auto high_res = baysor::louvain_partition(adj, /*resolution=*/2.0, /*max_passes=*/100);

    std::set<int> low_groups(low_res.begin(), low_res.end());
    std::set<int> high_groups(high_res.begin(), high_res.end());
    EXPECT_EQ(low_groups.size(), 3u);
    EXPECT_EQ(high_groups.size(), 6u);
}

TEST(MoleculeClustering, LouvainOnLocalGraphTracksLocalPatchesNotRepeatedLatentTypes) {
    // Six spatial patches arranged in a chain. Patches (0,3), (1,4), (2,5) are
    // intended to represent repeated latent types, but the graph only contains
    // local spatial edges, so Louvain can only recover local communities.
    constexpr int patches = 6;
    constexpr int patch_size = 5;
    constexpr int n = patches * patch_size;
    std::vector<int> src, dst;
    std::vector<double> wts;

    auto add_clique = [&](int start, int size, double w) {
        for (int i = start; i < start + size; ++i) {
            for (int j = i + 1; j < start + size; ++j) {
                src.push_back(i);
                dst.push_back(j);
                wts.push_back(w);
            }
        }
    };

    for (int p = 0; p < patches; ++p) {
        add_clique(p * patch_size, patch_size, 1.0);
    }
    for (int p = 0; p < patches - 1; ++p) {
        src.push_back(p * patch_size);
        dst.push_back((p + 1) * patch_size);
        wts.push_back(0.02);
    }

    auto adj = baysor::AdjList::from_edge_list(
        src.data(), dst.data(), wts.data(), static_cast<int>(src.size()), n);
    auto membership = baysor::louvain_partition(adj, /*resolution=*/1.0, /*max_passes=*/100);

    std::set<int> groups(membership.begin(), membership.end());
    EXPECT_EQ(groups.size(), 6u);
}

TEST(MoleculeClustering, KnnSimilarityGraphConnectsRepeatedLatentTypes) {
    constexpr int patches = 6;
    constexpr int patch_size = 5;
    constexpr int n = patches * patch_size;

    Eigen::MatrixXf mol_vecs(2, n);
    std::vector<double> confidence(n, 1.0);
    for (int p = 0; p < patches; ++p) {
        Eigen::Vector2f center;
        switch (p % 3) {
            case 0: center << 1.0f, 0.1f; break;
            case 1: center << 0.1f, 1.0f; break;
            default: center << -1.0f, 0.1f; break;
        }
        for (int i = 0; i < patch_size; ++i) {
            const int idx = p * patch_size + i;
            mol_vecs.col(idx) = center;
            mol_vecs(0, idx) += 0.01f * static_cast<float>(i);
            mol_vecs(1, idx) -= 0.01f * static_cast<float>(i);
        }
    }

    auto adj = baysor::build_knn_similarity_graph(mol_vecs, confidence, /*k=*/4);
    auto membership = baysor::louvain_partition(adj, /*resolution=*/1.0, /*max_passes=*/100);

    std::set<int> groups(membership.begin(), membership.end());
    EXPECT_EQ(groups.size(), 3u);
    for (int i = 0; i < patch_size; ++i) {
        EXPECT_EQ(membership[i], membership[3 * patch_size + i]);
        EXPECT_EQ(membership[patch_size + i], membership[4 * patch_size + i]);
        EXPECT_EQ(membership[2 * patch_size + i], membership[5 * patch_size + i]);
    }
}

TEST(MoleculeClustering, GraphPartitionToTargetReturnsRequestedCountForLouvainAndLeiden) {
    constexpr int patches = 6;
    constexpr int patch_size = 5;
    constexpr int n = patches * patch_size;

    Eigen::MatrixXf mol_vecs(2, n);
    std::vector<double> confidence(n, 1.0);
    for (int p = 0; p < patches; ++p) {
        Eigen::Vector2f center;
        switch (p % 3) {
            case 0: center << 1.0f, 0.1f; break;
            case 1: center << 0.1f, 1.0f; break;
            default: center << -1.0f, 0.1f; break;
        }
        for (int i = 0; i < patch_size; ++i) {
            const int idx = p * patch_size + i;
            mol_vecs.col(idx) = center;
            mol_vecs(0, idx) += 0.01f * static_cast<float>(i);
            mol_vecs(1, idx) -= 0.01f * static_cast<float>(i);
        }
    }

    auto adj = baysor::build_knn_similarity_graph(mol_vecs, confidence, /*k=*/4);

    baysor::GraphClusteringSummary lou_summary;
    auto louvain = baysor::graph_partition_to_target(
        adj, mol_vecs, confidence, baysor::ClusterMethod::Louvain,
        /*target_clusters=*/3, /*resolution_seed=*/1.0, /*max_passes=*/100, &lou_summary
    );
    baysor::GraphClusteringSummary lei_summary;
    auto leiden = baysor::graph_partition_to_target(
        adj, mol_vecs, confidence, baysor::ClusterMethod::Leiden,
        /*target_clusters=*/3, /*resolution_seed=*/1.0, /*max_passes=*/100, &lei_summary
    );

    std::set<int> lou_groups(louvain.begin(), louvain.end());
    std::set<int> lei_groups(leiden.begin(), leiden.end());
    EXPECT_EQ(lou_groups.size(), 3u);
    EXPECT_EQ(lei_groups.size(), 3u);
    EXPECT_EQ(lou_summary.final_clusters, 3);
    EXPECT_EQ(lei_summary.final_clusters, 3);
}

TEST(MoleculeClustering, GraphPartitionToTargetStableAcrossThreadCounts) {
    constexpr int patches = 8;
    constexpr int patch_size = 6;
    constexpr int n = patches * patch_size;

    Eigen::MatrixXf mol_vecs(3, n);
    std::vector<double> confidence(n, 1.0);
    for (int p = 0; p < patches; ++p) {
        Eigen::Vector3f center;
        switch (p % 4) {
            case 0: center << 1.0f, 0.1f, 0.0f; break;
            case 1: center << 0.1f, 1.0f, 0.1f; break;
            case 2: center << -1.0f, 0.2f, 0.0f; break;
            default: center << 0.2f, -1.0f, 0.1f; break;
        }
        for (int i = 0; i < patch_size; ++i) {
            const int idx = p * patch_size + i;
            mol_vecs.col(idx) = center;
            mol_vecs(0, idx) += 0.01f * static_cast<float>(i);
            mol_vecs(1, idx) -= 0.005f * static_cast<float>(i);
            confidence[idx] = 0.8 + 0.01 * static_cast<double>(i);
        }
    }

    auto adj = baysor::build_knn_similarity_graph(mol_vecs, confidence, /*k=*/5);

    for (auto method : {baysor::ClusterMethod::Louvain, baysor::ClusterMethod::Leiden}) {
        int old_threads = baysor::thread_pool_size();
        baysor::set_thread_pool_size(1);
        baysor::GraphClusteringSummary summary_1;
        auto assignment_1 = baysor::graph_partition_to_target(
            adj, mol_vecs, confidence, method,
            /*target_clusters=*/4, /*resolution_seed=*/1.0, /*max_passes=*/100, &summary_1
        );
        baysor::set_thread_pool_size(4);
        baysor::GraphClusteringSummary summary_4;
        auto assignment_4 = baysor::graph_partition_to_target(
            adj, mol_vecs, confidence, method,
            /*target_clusters=*/4, /*resolution_seed=*/1.0, /*max_passes=*/100, &summary_4
        );
        baysor::set_thread_pool_size(old_threads);

        EXPECT_EQ(assignment_4, assignment_1);
        EXPECT_EQ(summary_4.micro_clusters, summary_1.micro_clusters);
        EXPECT_EQ(summary_4.final_clusters, summary_1.final_clusters);
        EXPECT_DOUBLE_EQ(summary_4.chosen_resolution, summary_1.chosen_resolution);
        ASSERT_EQ(summary_4.move_fracs.size(), summary_1.move_fracs.size());
        for (size_t i = 0; i < summary_1.move_fracs.size(); ++i) {
            EXPECT_DOUBLE_EQ(summary_4.move_fracs[i], summary_1.move_fracs[i]);
        }
    }
}

TEST(MoleculeClustering, DispatcherReturnsEmptyForNoneMethod) {
    Eigen::MatrixXd pos(2, 2);
    pos <<
        0.0, 1.0,
        0.0, 0.0;
    std::vector<int> genes = {1, 2};
    std::vector<double> confidence = {1.0, 1.0};
    const int edge_src[] = {0};
    const int edge_dst[] = {1};
    const double edge_wt[] = {1.0};
    auto adj = baysor::AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 1, 2);

    baysor::ClusteringOptions opts;
    opts.method = baysor::ClusterMethod::None;

    auto result = baysor::cluster_molecules(pos, genes, adj, confidence, opts, /*verbose=*/false);
    EXPECT_TRUE(result.assignment.empty());
    EXPECT_TRUE(result.diffs.empty());
}

// ============================================================================
// Neighborhood composition
// ============================================================================

TEST(NeighborhoodComposition, CountMatrix) {
    // 6 molecules, 3 genes, on a line
    Eigen::MatrixXd pos(2, 6);
    for (int i = 0; i < 6; ++i) {
        pos(0, i) = i * 1.0;
        pos(1, i) = 0.0;
    }
    std::vector<int> genes = {1, 2, 3, 1, 2, 3};

    auto cm = baysor::neighborhood_count_matrix(pos, genes, 3, 3);

    // Should be 3 (genes) x 6 (molecules)
    EXPECT_EQ(cm.rows(), 3);
    EXPECT_EQ(cm.cols(), 6);

    // Each column should have non-zero entries (neighbors exist)
    for (int i = 0; i < 6; ++i) {
        double col_sum = 0;
        for (int j = 0; j < 3; ++j) {
            col_sum += cm.coeff(j, i);
        }
        EXPECT_GT(col_sum, 0.0);
    }
}

TEST(NeighborhoodComposition, EstimateGeneVectors) {
    // Create a small count matrix (3 genes x 10 molecules)
    Eigen::SparseMatrix<float> cm(3, 10);
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    for (int j = 0; j < 10; ++j) {
        for (int i = 0; i < 3; ++i) {
            if (dist(rng) > 0.3f) {
                cm.insert(i, j) = dist(rng);
            }
        }
    }
    cm.makeCompressed();

    std::vector<int> gene_ids = {1, 2, 3, 1, 2, 3, 1, 2, 3, 1};

    auto gene_vecs = baysor::estimate_gene_vectors(cm, gene_ids, 2, "ri", false);

    // Should be 2 (components) x 3 (genes)
    EXPECT_EQ(gene_vecs.rows(), 2);
    EXPECT_EQ(gene_vecs.cols(), 3);
}

TEST(NeighborhoodComposition, EstimateGeneVectorsPerMolecule) {
    Eigen::SparseMatrix<float> cm(3, 10);
    for (int j = 0; j < 10; ++j) {
        cm.insert(j % 3, j) = 1.0f;
    }
    cm.makeCompressed();

    std::vector<int> gene_ids(10, 1);
    auto mol_vecs = baysor::estimate_gene_vectors(cm, gene_ids, 5, "ri", true);

    // Should be 5 (components) x 10 (molecules)
    EXPECT_EQ(mol_vecs.rows(), 5);
    EXPECT_EQ(mol_vecs.cols(), 10);
}

TEST(NeighborhoodComposition, ProjectNeighborhoodVectorsMatchesLoggedCountMatrixProjection) {
    Eigen::MatrixXd pos(2, 8);
    for (int i = 0; i < 8; ++i) {
        pos(0, i) = static_cast<double>(i);
        pos(1, i) = static_cast<double>(i % 2);
    }
    std::vector<int> genes = {1, 2, 3, 1, 2, 3, 1, 2};
    std::vector<int> query_ids = {1, 3, 4, 6};

    const int k_neighbors = 4;
    const int n_genes = 3;
    const double dist_floor = baysor::neighborhood_distance_floor(pos, &query_ids);

    auto cm = baysor::neighborhood_count_matrix_subset(
        pos, genes, query_ids, k_neighbors, n_genes, nullptr, true, true, dist_floor
    );
    for (int k = 0; k < cm.outerSize(); ++k) {
        for (Eigen::SparseMatrix<float>::InnerIterator it(cm, k); it; ++it) {
            it.valueRef() = static_cast<float>(std::log(it.value() * 10000.0f + 1e-5f));
        }
    }

    auto gene_emb_t = baysor::estimate_gene_vectors(cm, genes, 4, "ri", false);
    auto expected = baysor::project_gene_vectors(gene_emb_t, cm);
    auto actual = baysor::project_neighborhood_vectors(
        pos, genes, k_neighbors, gene_emb_t, n_genes, &query_ids, nullptr, true, true, dist_floor, true
    );

    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());
    EXPECT_LT((actual - expected).norm(), 1e-4f);
}

// ============================================================================
// Color embedding
// ============================================================================

TEST(ColorUtils, NormalizeEmbedding) {
    Eigen::MatrixXd emb(3, 100);
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 100; ++j) {
            emb(i, j) = dist(rng);
        }
    }

    baysor::normalize_embedding_to_lab_range(emb);

    // L channel should be in [~10, ~90]
    for (int j = 0; j < 100; ++j) {
        EXPECT_GE(emb(0, j), 0.0);
        EXPECT_LE(emb(0, j), 100.0);
    }
    // a, b channels should be in [-100, 100]
    for (int j = 0; j < 100; ++j) {
        EXPECT_GE(emb(1, j), -110.0);
        EXPECT_LE(emb(1, j), 110.0);
        EXPECT_GE(emb(2, j), -110.0);
        EXPECT_LE(emb(2, j), 110.0);
    }
}

TEST(ColorUtils, EmbeddingToHex) {
    // Middle gray in LAB: L=50, a=0, b=0
    Eigen::MatrixXd lab(3, 2);
    lab.col(0) << 50.0, 0.0, 0.0;  // neutral gray
    lab.col(1) << 90.0, 0.0, 0.0;  // light gray

    auto hexes = baysor::embedding_to_hex(lab);

    ASSERT_EQ(hexes.size(), 2u);
    // Should start with '#' and be 7 chars
    EXPECT_EQ(hexes[0].size(), 7u);
    EXPECT_EQ(hexes[0][0], '#');
    EXPECT_EQ(hexes[1].size(), 7u);
    EXPECT_EQ(hexes[1][0], '#');
}

TEST(ColorUtils, GeneCompositionColorEmbedding) {
    // Create synthetic mol_vecs: 10 components x 50 molecules
    Eigen::MatrixXf mol_vecs(10, 50);
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 50; ++j) {
            mol_vecs(i, j) = static_cast<float>(dist(rng));
        }
    }
    // The embedding helper fits on high-confidence molecules (>= 0.95).
    std::vector<double> confidence(50, 0.99);

    auto colors = baysor::gene_composition_color_embedding(mol_vecs, confidence, 30);

    ASSERT_EQ(colors.size(), 50u);
    for (const auto& c : colors) {
        EXPECT_EQ(c.size(), 7u);
        EXPECT_EQ(c[0], '#');
    }
}

TEST(ColorUtils, GeneCompositionColorEmbeddingFallsBackWithoutHighConfidenceAnchors) {
    Eigen::MatrixXf mol_vecs(10, 8);
    std::mt19937 rng(7);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < mol_vecs.rows(); ++i) {
        for (int j = 0; j < mol_vecs.cols(); ++j) {
            mol_vecs(i, j) = static_cast<float>(dist(rng));
        }
    }

    std::vector<double> confidence(mol_vecs.cols(), 0.20);
    auto colors = baysor::gene_composition_color_embedding(mol_vecs, confidence, 30);

    ASSERT_EQ(colors.size(), static_cast<size_t>(mol_vecs.cols()));
    for (const auto& c : colors) {
        EXPECT_EQ(c, "#808080");
    }
}

TEST(ColorUtils, GeneCompositionColorEmbeddingLowersThresholdWhenNeeded) {
    Eigen::MatrixXf mol_vecs(10, 50);
    std::mt19937 rng(11);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < mol_vecs.rows(); ++i) {
        for (int j = 0; j < mol_vecs.cols(); ++j) {
            mol_vecs(i, j) = static_cast<float>(dist(rng));
        }
    }

    // No molecule reaches the historical 0.95 cutoff, but all molecules qualify
    // once the adaptive ladder relaxes to 0.90.
    std::vector<double> confidence(mol_vecs.cols(), 0.92);
    auto colors = baysor::gene_composition_color_embedding(mol_vecs, confidence, 30);

    ASSERT_EQ(colors.size(), static_cast<size_t>(mol_vecs.cols()));
    bool saw_non_fallback = false;
    for (const auto& c : colors) {
        EXPECT_EQ(c.size(), 7u);
        EXPECT_EQ(c[0], '#');
        if (c != "#808080") saw_non_fallback = true;
    }
    EXPECT_TRUE(saw_non_fallback);
}

TEST(ColorUtils, GeneCompositionColorEmbeddingStreamingDoesNotCollapseToSingleColor) {
    constexpr int groups = 3;
    constexpr int per_group = 30;
    constexpr int n = groups * per_group;

    Eigen::MatrixXd pos(2, n);
    std::vector<int> genes(n, 1);
    std::vector<double> confidence(n, 0.99);

    for (int g = 0; g < groups; ++g) {
        for (int i = 0; i < per_group; ++i) {
            const int idx = g * per_group + i;
            pos(0, idx) = 50.0 * g + 0.25 * i;
            pos(1, idx) = 10.0 * (i % 5);
            genes[idx] = 1 + g;
        }
    }

    auto colors = baysor::gene_composition_color_embedding_streaming(
        pos, genes, /*n_genes=*/3, confidence,
        /*k_neighbors=*/8,
        /*basis_sample_size=*/60,
        /*sample_size=*/30,
        /*seed=*/1,
        /*n_pca_dims=*/3
    );

    ASSERT_EQ(colors.size(), static_cast<size_t>(n));
    std::set<std::string> unique(colors.begin(), colors.end());
    EXPECT_GT(unique.size(), 1u);
    for (const auto& c : colors) {
        EXPECT_EQ(c.size(), 7u);
        EXPECT_EQ(c[0], '#');
    }
}

TEST(ColorUtils, GeneCompositionColorEmbeddingStreamingMatchesPrecomputedModel) {
    constexpr int groups = 3;
    constexpr int per_group = 24;
    constexpr int n = groups * per_group;

    Eigen::MatrixXd pos(2, n);
    std::vector<int> genes(n, 1);
    std::vector<double> confidence(n, 0.99);

    for (int g = 0; g < groups; ++g) {
        for (int i = 0; i < per_group; ++i) {
            const int idx = g * per_group + i;
            pos(0, idx) = 30.0 * g + 0.5 * i;
            pos(1, idx) = 5.0 * (i % 6);
            genes[idx] = 1 + g;
        }
    }

    auto expected = baysor::gene_composition_color_embedding_streaming(
        pos, genes, /*n_genes=*/3, confidence,
        /*k_neighbors=*/8,
        /*basis_sample_size=*/48,
        /*sample_size=*/24,
        /*seed=*/7,
        /*n_pca_dims=*/3
    );

    auto model = baysor::fit_ncv_projected_model(
        pos, genes, /*n_genes=*/3, confidence,
        /*k_neighbors=*/8,
        /*basis_sample_size=*/48,
        /*n_components=*/20,
        /*include_full_projection=*/true
    );

    auto actual = baysor::gene_composition_color_embedding_streaming(
        pos, genes, /*n_genes=*/3, confidence,
        /*k_neighbors=*/8,
        /*basis_sample_size=*/48,
        /*sample_size=*/24,
        /*seed=*/7,
        /*n_pca_dims=*/3,
        /*graph_k=*/15,
        &model
    );

    EXPECT_EQ(actual, expected);
}

// ============================================================================
// Main BMM loop
// ============================================================================

TEST(BmmLoop, StableDisconnectedComponentsRemainStable) {
    auto data = make_disconnected_bmm_data();

    baysor::bmm(data,
                /*min_molecules_drop=*/2,
                /*n_iters=*/3,
                /*assignment_history_depth=*/3,
                /*verbose=*/false,
                /*component_split_step=*/3,
                /*refine=*/false,
                /*freeze_composition=*/false,
                /*freeze_position=*/false,
                /*freeze_components=*/false,
                /*tol=*/0.0,
                /*min_molecules_display=*/2);

    EXPECT_EQ(data.n_components(), 3);
    EXPECT_EQ(data.assignment, (std::vector<int>{1, 1, 1, 2, 2, 2, 3, 3}));
    EXPECT_EQ(data.assignment_history.size(), 3u);
    EXPECT_EQ(data.n_components_trace.size(), 4u);  // initial trace + one per iteration
}

TEST(BmmLoop, PriorSegmentationAdjustmentCanFlipWinner) {
    auto no_prior = make_two_component_competition_data();
    auto with_prior = make_two_component_competition_data();

    with_prior.segment_per_molecule = {1, 1, 1, 1, 1, 2};
    with_prior.n_molecules_per_segment = {5, 1};
    with_prior.prior_seg_confidence = 0.5;
    with_prior.update_n_mols_per_segment();

    baysor::expect_dirichlet_spatial(no_prior, /*stochastic=*/false);
    baysor::expect_dirichlet_spatial(with_prior, /*stochastic=*/false);

    EXPECT_EQ(no_prior.assignment[0], 2);
    EXPECT_EQ(with_prior.assignment[0], 1);
}

TEST(BmmLoop, EstepStatsReportsChangedAssignments) {
    auto data = make_two_component_competition_data();
    auto before = data.assignment;

    auto stats = baysor::expect_dirichlet_spatial(data, /*stochastic=*/false);

    std::int64_t expected_changed = 0;
    for (size_t i = 0; i < before.size(); ++i) {
        if (before[i] != data.assignment[i]) ++expected_changed;
    }
    EXPECT_EQ(stats.n_changed, expected_changed);
}

TEST(BmmLoop, ConnectedComponentSplitMatchesAcrossThreadCounts) {
    auto one_thread = make_disconnected_bmm_data();
    auto many_threads = make_disconnected_bmm_data();
    one_thread.assignment.assign(one_thread.assignment.size(), 1);
    many_threads.assignment = one_thread.assignment;

    int old_threads = baysor::thread_pool_size();
    baysor::set_thread_pool_size(1);
    baysor::split_cells_by_connected_components(one_thread);
    baysor::set_thread_pool_size(4);
    baysor::split_cells_by_connected_components(many_threads);
    baysor::set_thread_pool_size(old_threads);

    const std::vector<int> expected = {1, 1, 1, 0, 0, 0, 0, 0};
    EXPECT_EQ(one_thread.assignment, expected);
    EXPECT_EQ(many_threads.assignment, expected);
}

TEST(BmmLoop, ClusterPenaltyCanFlipWinner) {
    auto no_penalty = make_two_component_competition_data();
    auto with_penalty = make_two_component_competition_data();

    no_penalty.cluster_per_cell = {1, 2, 3};
    no_penalty.cluster_per_molecule = {1, 1, 1, 1, 2, 3};
    no_penalty.cluster_penalty_mult = 1.0;

    with_penalty.cluster_per_cell = {1, 2, 3};
    with_penalty.cluster_per_molecule = {1, 1, 1, 1, 2, 3};
    with_penalty.cluster_penalty_mult = 0.25;

    baysor::expect_dirichlet_spatial(no_penalty, /*stochastic=*/false);
    baysor::expect_dirichlet_spatial(with_penalty, /*stochastic=*/false);

    EXPECT_EQ(no_penalty.assignment[0], 2);
    EXPECT_EQ(with_penalty.assignment[0], 1);
}

TEST(BmmLoop, ClusterPenaltySkipsLastComponentLikeJulia) {
    auto baseline = make_last_component_competition_data();
    auto julia_like = make_last_component_competition_data();

    julia_like.cluster_per_cell = {1, 2, 3};
    julia_like.cluster_per_molecule = {1, 1, 1, 1, 3, 2};
    julia_like.cluster_penalty_mult = 0.25;

    baysor::expect_dirichlet_spatial(baseline, /*stochastic=*/false);
    baysor::expect_dirichlet_spatial(julia_like, /*stochastic=*/false);

    EXPECT_EQ(baseline.assignment[0], 3);
    EXPECT_EQ(julia_like.assignment[0], 3);
}

TEST(BmmLoop, NoiseDensityCompetesWithComponents) {
    auto low_noise = make_noise_competition_data();
    auto high_noise = make_noise_competition_data();

    low_noise.noise_density = 1e-6;
    high_noise.noise_density = 1.0;

    baysor::expect_dirichlet_spatial(low_noise, /*stochastic=*/false);
    baysor::expect_dirichlet_spatial(high_noise, /*stochastic=*/false);

    EXPECT_EQ(low_noise.assignment[0], 1);
    EXPECT_EQ(high_noise.assignment[0], 0);
}

TEST(BmmLoop, NoiseDensitySkipsEmptyComponentsLikeJulia) {
    auto data = make_noise_density_skip_empty_component_data();

    baysor::maximize(data, /*freeze_composition=*/false, /*freeze_position=*/true);

    EXPECT_NEAR(data.noise_density, 1.0, 1e-12);
    EXPECT_EQ(data.components[0].composition_params.n_genes, 2);
    EXPECT_LE(data.components[1].composition_params.sum_counts, 1e-12);
}

TEST(BmmLoop, HigherDropThresholdPrunesStableSmallComponent) {
    auto julia_like = make_disconnected_bmm_data();
    auto cpp_like = make_disconnected_bmm_data();

    baysor::bmm(julia_like,
                /*min_molecules_drop=*/2,
                /*n_iters=*/1,
                /*assignment_history_depth=*/0,
                /*verbose=*/false,
                /*component_split_step=*/3,
                /*refine=*/false,
                /*freeze_composition=*/false,
                /*freeze_position=*/false,
                /*freeze_components=*/false,
                /*tol=*/0.0,
                /*min_molecules_display=*/2);

    baysor::bmm(cpp_like,
                /*min_molecules_drop=*/3,
                /*n_iters=*/1,
                /*assignment_history_depth=*/0,
                /*verbose=*/false,
                /*component_split_step=*/3,
                /*refine=*/false,
                /*freeze_composition=*/false,
                /*freeze_position=*/false,
                /*freeze_components=*/false,
                /*tol=*/0.0,
                /*min_molecules_display=*/3);

    EXPECT_EQ(julia_like.n_components(), 3);
    EXPECT_EQ(julia_like.assignment, (std::vector<int>{1, 1, 1, 2, 2, 2, 3, 3}));

    EXPECT_EQ(cpp_like.n_components(), 2);
    EXPECT_EQ(cpp_like.assignment, (std::vector<int>{1, 1, 1, 2, 2, 2, 0, 0}));
}

TEST(BmmLoop, PositiveToleranceStopsBeforeConfiguredIterations) {
    auto fixed_iters = make_disconnected_bmm_data();
    auto early_stop = make_disconnected_bmm_data();

    constexpr int n_iters = 30;

    baysor::bmm(fixed_iters,
                /*min_molecules_drop=*/2,
                /*n_iters=*/n_iters,
                /*assignment_history_depth=*/n_iters,
                /*verbose=*/false,
                /*component_split_step=*/3,
                /*refine=*/false,
                /*freeze_composition=*/false,
                /*freeze_position=*/false,
                /*freeze_components=*/false,
                /*tol=*/0.0,
                /*min_molecules_display=*/2);

    baysor::bmm(early_stop,
                /*min_molecules_drop=*/2,
                /*n_iters=*/n_iters,
                /*assignment_history_depth=*/n_iters,
                /*verbose=*/false,
                /*component_split_step=*/3,
                /*refine=*/false,
                /*freeze_composition=*/false,
                /*freeze_position=*/false,
                /*freeze_components=*/false,
                /*tol=*/0.01,
                /*min_molecules_display=*/2);

    EXPECT_EQ(fixed_iters.assignment_history.size(), static_cast<size_t>(n_iters));
    EXPECT_EQ(fixed_iters.n_components_trace.size(), static_cast<size_t>(n_iters + 1));

    EXPECT_LT(early_stop.assignment_history.size(), static_cast<size_t>(n_iters));
    EXPECT_LT(early_stop.n_components_trace.size(), static_cast<size_t>(n_iters + 1));
    EXPECT_EQ(early_stop.assignment, fixed_iters.assignment);
}

TEST(RunReport, GeneratesDiagnosticHtml) {
    baysor::MoleculeData data;
    data.x = {0.0, 1.0, 2.0, 3.0};
    data.y = {0.0, 0.0, 1.0, 1.0};
    data.gene = {1, 1, 1, 1};
    data.gene_names = {"GeneA"};
    data.confidence = {0.95, 0.85, 0.25, 0.10};
    data.prior_segmentation = {1, 1, 0, 2};

    std::vector<double> edge_lengths = {0.1, 0.12, 1.0, 1.2};
    baysor::NoiseFitResult noise_result;
    noise_result.assignment_probs = Eigen::MatrixXd::Zero(4, 2);
    noise_result.assignment = {1, 1, 2, 2};
    noise_result.signal_mu = 0.11;
    noise_result.signal_sigma = 0.02;
    noise_result.noise_mu = 1.1;
    noise_result.noise_sigma = 0.1;

    std::vector<int> assignment = {1, 1, 0, 2};
    std::vector<std::unordered_map<int, int>> trace = {
        {{1, 2}, {2, 1}},
        {{1, 2}, {2, 2}},
        {{1, 2}, {2, 2}}
    };
    std::vector<double> assignment_conf = {1.0, 1.0, 0.0, 1.0};

    Eigen::MatrixXd cell_stats(2, 5);
    cell_stats <<
        0.5, 0.0, 2.0, 2.0, 1.1,
        3.0, 1.0, 2.5, 1.2, 1.3;
    std::vector<std::string> cell_cols = {"x", "y", "area", "density", "elongation"};

    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Column;
    prior.column_name = "cell_id";

    baysor::ClusteringResult clustering_result;
    clustering_result.assignment = {1, 1, 2, 2};

    baysor::NcvReportEmbedding ncv_report;
    ncv_report.colors = {"#ff0000", "#00ff00", "#0000ff", "#ff00ff"};
    ncv_report.sample_ids = {0, 2, 3};
    ncv_report.sample_umap_x = {0.0, 1.0, 2.0};
    ncv_report.sample_umap_y = {0.0, 1.0, 0.5};

    auto html = baysor::generate_run_diagnostic_html(
        data, edge_lengths, noise_result, 10, assignment, trace, assignment_conf,
        &clustering_result, &ncv_report, cell_stats, cell_cols, prior, 4.5, "1.2");

    EXPECT_NE(html.find("Baysor Run Report"), std::string::npos);
    EXPECT_NE(html.find("Segmentation convergence"), std::string::npos);
    EXPECT_NE(html.find("Molecule confidence"), std::string::npos);
    EXPECT_NE(html.find("Prior column: cell_id"), std::string::npos);
    EXPECT_NE(html.find("NCV / clustering manifold"), std::string::npos);
    EXPECT_NE(html.find("Colored by cluster assignment"), std::string::npos);
}

TEST(RunReport, GeneratesSegmentationHtml) {
    baysor::MoleculeData data;
    data.x = {0.0, 1.0, 2.0, 3.0};
    data.y = {0.0, 0.0, 1.0, 1.0};
    data.gene = {1, 1, 1, 1};
    data.gene_names = {"GeneA"};

    std::vector<int> assignment = {1, 1, 0, 2};
    std::vector<std::string> ncv_color = {"#ff0000", "#00ff00", "#0000ff", "#ff00ff"};
    std::vector<int> mol_clusters = {1, 1, 2, 2};

    baysor::PolygonCollection polygons;
    Eigen::MatrixXd poly_2xN(2, 4);
    poly_2xN << 0.0, 2.0, 2.0, 0.0,
                0.0, 0.0, 1.5, 1.5;
    polygons["cell_1"] = poly_2xN;

    Eigen::MatrixXd poly_Nx2(4, 2);
    poly_Nx2 << 0.5, 0.2,
                1.0, 0.2,
                1.0, 0.8,
                0.5, 0.8;
    polygons["cell_2"] = poly_Nx2;

    auto html = baysor::generate_run_segmentation_html(
        data, assignment, ncv_color, &mol_clusters, &polygons);

    EXPECT_NE(html.find("Final cell assignment"), std::string::npos);
    EXPECT_NE(html.find("Local expression similarity (NCV)"), std::string::npos);
    EXPECT_NE(html.find("Molecule clustering"), std::string::npos);
    EXPECT_NE(html.find("data:image/png;base64,"), std::string::npos);
}

TEST(NeighborhoodComposition, EstimateGeneVectorsMatchesReferencePerMolecule) {
    Eigen::SparseMatrix<float> count_matrix = make_random_count_matrix(64, 320, 18, 123);
    std::vector<int> gene_ids(count_matrix.cols(), 1);

    Eigen::MatrixXf actual = baysor::estimate_gene_vectors(count_matrix, gene_ids, 8, "ri", true);
    Eigen::MatrixXf expected = estimate_gene_vectors_reference(count_matrix, 8, true);

    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());

    float denom = std::max(1e-6f, expected.norm());
    float rel_diff = (actual - expected).norm() / denom;
    EXPECT_LT(rel_diff, 1e-4f);
}

TEST(NeighborhoodComposition, EstimateGeneVectorsMatchesReferencePerGene) {
    Eigen::SparseMatrix<float> count_matrix = make_random_count_matrix(48, 180, 12, 456);
    std::vector<int> gene_ids(count_matrix.cols(), 1);

    Eigen::MatrixXf actual = baysor::estimate_gene_vectors(count_matrix, gene_ids, 6, "ri", false);
    Eigen::MatrixXf expected = estimate_gene_vectors_reference(count_matrix, 6, false);

    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());

    float denom = std::max(1e-6f, expected.norm());
    float rel_diff = (actual - expected).norm() / denom;
    EXPECT_LT(rel_diff, 1e-4f);
}
