// BUG-6 (review must-fix): integer config keys must accept integral floats
// and TOML digit separators, like Julia.
//
// Julia (Baysor v0.7.1, src/utils/options.jl): integer fields such as
// `min_molecules_per_cell::Int` are filled from the parsed TOML by
// Configurations.jl, which converts with `Base.convert` -- so
// `min_molecules_per_cell = 50.0` becomes 50 (convert(Int, 50.0) == 50) and
// `1e2` becomes 100. TOML.jl follows the TOML number grammar, so `1_000`
// parses as 1000 (digit separators) and as a double in float keys. The BUG-3
// strictness went too far and rejected all of these with "expected an
// integer", breaking configs Julia (and the old C++) accepted.
//
// Still errors, like Julia's InexactError / TOML parse errors: non-integral
// values (50.5), values out of `int` range, and garbage. Every error must
// name the key, the value and the expected type (see Bug3_ConfigErrors in
// tests/test_bugfix_correctness.cpp for the garbage cases).

#include <gtest/gtest.h>

#include "baysor/utils/options.h"

#include <fstream>
#include <stdexcept>
#include <string>

#include "test_cov_helpers.h"

namespace {

using baysor_test::TempDir;

// Load `content` as a TOML config; propagates load_config's runtime_error.
baysor::RunOptions load_config_text(const std::string& content) {
    TempDir dir("bug6_cfg");
    const std::string p = dir.file("cfg.toml");
    {
        std::ofstream f(p);
        f << content;
    }
    return baysor::load_config(p);
}

// Load `content` as a TOML config; return the runtime_error message, or an
// empty string if load_config did not throw.
std::string config_error(const std::string& content) {
    try {
        (void)load_config_text(content);
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return std::string();
}

}  // namespace

// ---------------------------------------------------------------------------
// Accepted (Julia converts / TOML parses these as integers)
// ---------------------------------------------------------------------------

TEST(Bug6_ConfigIntValues, IntegralFloat50IsAcceptedAsInt50) {
    // Configurations.jl: Base.convert(Int, 50.0) == 50.
    const auto opts = load_config_text(
        "[molecules]\nmin_molecules_per_cell = 50.0\n");
    EXPECT_EQ(opts.molecules.min_molecules_per_cell, 50);
}

TEST(Bug6_ConfigIntValues, ExponentForm1e2IsAcceptedAsInt100) {
    const auto opts = load_config_text(
        "[molecules]\nmin_molecules_per_cell = 1e2\n");
    EXPECT_EQ(opts.molecules.min_molecules_per_cell, 100);
}

TEST(Bug6_ConfigIntValues, DigitSeparator1000IsAcceptedAsInt1000) {
    // TOML.jl: `1_000` is 1000 (underscores between digits are separators).
    const auto opts = load_config_text(
        "[molecules]\nmin_molecules_per_cell = 1_000\n");
    EXPECT_EQ(opts.molecules.min_molecules_per_cell, 1000);
}

TEST(Bug6_ConfigIntValues, DigitSeparatorInDoubleKeyIsAccepted) {
    // Same TOML grammar applies to float keys.
    const auto opts = load_config_text(
        "[molecules]\nx_min = 1_234.5\n");
    EXPECT_DOUBLE_EQ(opts.molecules.x_min, 1234.5);
}

// ---------------------------------------------------------------------------
// Still rejected (Julia: InexactError / TOML parse error)
// ---------------------------------------------------------------------------

TEST(Bug6_ConfigIntValues, NonIntegralFloatIsRejected) {
    const std::string msg = config_error(
        "[molecules]\nmin_molecules_per_cell = 50.5\n");
    ASSERT_FALSE(msg.empty()) << "non-integral value silently kept default";
    EXPECT_NE(msg.find("min_molecules_per_cell"), std::string::npos) << msg;
    EXPECT_NE(msg.find("50.5"), std::string::npos) << msg;
    EXPECT_NE(msg.find("integer"), std::string::npos) << msg;
}

TEST(Bug6_ConfigIntValues, IntegralFloatPastIntMaxIsRejected) {
    // A whole number, but one past std::numeric_limits<int>::max();
    // convert(Int, 2147483648.0) is an InexactError for an Int32 field.
    const std::string msg = config_error(
        "[segmentation]\nn_clusters = 2147483648.0\n");
    ASSERT_FALSE(msg.empty()) << "out-of-int-range value silently kept default";
    EXPECT_NE(msg.find("n_clusters"), std::string::npos) << msg;
    EXPECT_NE(msg.find("2147483648.0"), std::string::npos) << msg;
    EXPECT_NE(msg.find("integer"), std::string::npos) << msg;
}

TEST(Bug6_ConfigIntValues, MalformedSeparatorIsRejected) {
    // TOML requires separators to sit between digits; `1__000` is a parse
    // error in TOML.jl, so it must not be silently repaired here either.
    const std::string msg = config_error(
        "[molecules]\nmin_molecules_per_cell = 1__000\n");
    ASSERT_FALSE(msg.empty()) << "malformed separator silently kept default";
    EXPECT_NE(msg.find("min_molecules_per_cell"), std::string::npos) << msg;
    EXPECT_NE(msg.find("1__000"), std::string::npos) << msg;
    EXPECT_NE(msg.find("integer"), std::string::npos) << msg;
}
