// Delta-encoded assignment history (AssignmentHistory) and its readers:
// row reconstruction, trimming, component lifespans and the history vote
// must match row-based reference implementations exactly.
#include <gtest/gtest.h>

#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/processing/bmm_algorithm/tracing.h"
#include "baysor/processing/models/assignment_history.h"
#include "baysor/utils/thread_pool.h"

#include <Eigen/Dense>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <vector>

namespace {

using baysor::AssignmentHistory;
using baysor::BmmData;

class PoolSizeGuard {
public:
    explicit PoolSizeGuard(int n) : old_(baysor::thread_pool_size()) {
        baysor::set_thread_pool_size(n);
    }
    ~PoolSizeGuard() { baysor::set_thread_pool_size(old_); }

private:
    int old_;
};

// Rows with ~15 % of the molecules changing per step, guids in [0, n_guids).
std::vector<std::vector<int>> random_rows(int n_rows, int n, int n_guids, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> guid(0, n_guids - 1);
    std::bernoulli_distribution change(0.15);
    std::vector<std::vector<int>> rows(n_rows, std::vector<int>(n));
    for (int i = 0; i < n; ++i) rows[0][i] = guid(rng);
    for (int t = 1; t < n_rows; ++t) {
        for (int i = 0; i < n; ++i) rows[t][i] = change(rng) ? guid(rng) : rows[t - 1][i];
    }
    return rows;
}

// Row-based reference implementations of the delta-encoded readers.
std::pair<std::vector<int>, std::vector<double>> reference_vote(
    const std::vector<std::vector<int>>& rows, const std::vector<int>& component_guids
) {
    const int n = static_cast<int>(rows.back().size());
    std::unordered_map<int, int> guid_map;
    for (int i = 0; i < static_cast<int>(component_guids.size()); ++i) guid_map[component_guids[i]] = i + 1;

    std::vector<int> reassignment(n, 0);
    std::vector<double> match_frac(n, 0.0);
    for (int i = 0; i < n; ++i) {
        std::unordered_map<int, int> freq;
        int valid = 0;
        for (const auto& row : rows) {
            int g = row[i];
            if (g == 0 || guid_map.count(g)) {
                freq[g]++;
                valid++;
            }
        }
        int best_guid = 0, best_cnt = 0;
        for (auto& [g, c] : freq) {
            if (c > best_cnt) { best_cnt = c; best_guid = g; }
        }
        reassignment[i] = guid_map.count(best_guid) ? guid_map[best_guid] : 0;
        match_frac[i] = (valid > 0) ? static_cast<double>(best_cnt) / valid : 0.0;
    }
    return {reassignment, match_frac};
}

// Lifespan of a guid of the newest entry: the number of newest entries in a
// row that contain it.
std::unordered_map<int, int> reference_lifespan(const std::vector<std::vector<int>>& rows) {
    std::unordered_map<int, int> life;
    for (int g : rows.back()) {
        if (g == 0 || life.count(g)) continue;
        int t = static_cast<int>(rows.size()) - 1;
        while (t >= 0 && std::count(rows[t].begin(), rows[t].end(), g) > 0) --t;
        life[g] = static_cast<int>(rows.size()) - 1 - t;
    }
    return life;
}

BmmData<2> data_with_components(int n, const std::vector<int>& guids) {
    BmmData<2> data;
    data.position_data = Eigen::MatrixXd::Zero(2, n);
    data.composition_data.assign(n, 0);
    data.confidence.assign(n, 1.0);
    data.assignment.assign(n, 0);
    for (int g : guids) {
        baysor::MvNormal<2> pos(Eigen::Vector2d::Zero(), Eigen::Matrix2d::Identity());
        baysor::CategoricalSmoothed comp(2, 1.0);
        data.components.emplace_back(pos, comp, std::nullopt, g);
    }
    return data;
}

} // namespace

TEST(AssignmentHistory, RoundTripsRowsAndTrimsOldestEntries) {
    for (int n_threads : {1, 4}) {
        PoolSizeGuard guard(n_threads);
        const auto rows = random_rows(12, 50000, 9, 7u);
        AssignmentHistory h(rows);
        ASSERT_EQ(h.size(), rows.size());
        EXPECT_EQ(h.rows(), rows);
        EXPECT_EQ(h.back(), rows.back());

        // Deltas hold exactly the changed molecules, ascending, older values.
        for (size_t t = 0; t + 1 < rows.size(); ++t) {
            const auto& d = h.changes_after(t);
            size_t k = 0;
            for (int i = 0; i < 50000; ++i) {
                if (rows[t][i] == rows[t + 1][i]) continue;
                ASSERT_LT(k, d.size());
                EXPECT_EQ(d[k].mol, i);
                EXPECT_EQ(d[k].value, rows[t][i]);
                ++k;
            }
            EXPECT_EQ(k, d.size());
        }

        h.pop_front();
        h.pop_front();
        std::vector<std::vector<int>> rest(rows.begin() + 2, rows.end());
        EXPECT_EQ(h.rows(), rest);
        // Pushing after a pop reuses the evicted buffer and stays exact.
        h.push_back(rows[0]);
        rest.push_back(rows[0]);
        EXPECT_EQ(h.rows(), rest);

        while (!h.empty()) h.pop_front();
        EXPECT_TRUE(h.rows().empty());
    }
}

TEST(AssignmentHistory, RejectsRowsOfDifferentLength) {
    AssignmentHistory h;
    h.push_back({1, 2, 3});
    EXPECT_THROW(h.push_back({1, 2}), std::invalid_argument);
}

TEST(AssignmentHistory, LifespanMatchesRowBasedComputation) {
    for (unsigned seed : {1u, 2u, 3u}) {
        // Few guids that come and go, so streaks break at different depths.
        auto rows = random_rows(20, 300, 30, seed);
        EXPECT_EQ(baysor::estimate_component_lifespan(AssignmentHistory(rows)),
                  reference_lifespan(rows)) << "seed " << seed;
    }
    EXPECT_TRUE(baysor::estimate_component_lifespan(AssignmentHistory()).empty());
}

TEST(AssignmentHistory, HistoryVoteMatchesRowBasedComputation) {
    // Guid 0..6 occur in the history; guids 2 and 5 are not current, so they
    // are ignored; many ties exercise the hash-map iteration-order tie-break.
    const std::vector<int> current = {1, 3, 4, 6};
    for (int n_threads : {1, 4}) {
        PoolSizeGuard guard(n_threads);
        for (unsigned seed : {11u, 12u}) {
            const int n = 20000;
            auto rows = random_rows(10, n, 7, seed);
            auto data = data_with_components(n, current);
            data.assignment_history = AssignmentHistory(rows);

            auto [reassign, frac] = baysor::estimate_assignment_by_history(data);
            auto [ref_reassign, ref_frac] = reference_vote(rows, current);
            EXPECT_EQ(reassign, ref_reassign) << "threads " << n_threads << " seed " << seed;
            EXPECT_EQ(frac, ref_frac) << "threads " << n_threads << " seed " << seed;
        }
    }
}
