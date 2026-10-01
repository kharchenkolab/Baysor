// BUG-6: integer config keys accept integral floats and TOML digit
// separators, like Julia (Baysor v0.7.1): Configurations.jl converts with
// `Base.convert(Int, 50.0) == 50` and TOML.jl reads `1_000` as 1000.
// Non-integral values, values out of `int` range and malformed separators
// still fail (Julia: InexactError / TOML parse error), naming the key, the
// value and the expected type. Garbage values are covered by
// Bug3_ConfigErrors in tests/test_bugfix_correctness.cpp.

#include <gtest/gtest.h>

#include "baysor/utils/options.h"

#include <stdexcept>
#include <string>

#include "test_cov_helpers.h"

namespace {

baysor::RunOptions load_config_text(const std::string& content) {
    baysor_test::TempDir dir("bug6_cfg");
    return baysor::load_config(dir.write("cfg.toml", content));
}

}  // namespace

TEST(Bug6_ConfigIntValues, IntegralFloatsAndDigitSeparatorsAreAccepted) {
    const std::pair<const char*, int> int_cases[] = {
        {"50.0", 50}, {"1e2", 100}, {"1_000", 1000},
    };
    for (const auto& [value, expected] : int_cases) {
        const auto opts = load_config_text(std::string("[molecules]\nmin_molecules_per_cell = ") + value + "\n");
        EXPECT_EQ(opts.molecules.min_molecules_per_cell, expected) << value;
    }
    // The separator grammar applies to float keys as well.
    EXPECT_DOUBLE_EQ(load_config_text("[molecules]\nx_min = 1_234.5\n").molecules.x_min, 1234.5);
}

TEST(Bug6_ConfigIntValues, NonIntegralOutOfRangeAndMalformedValuesAreRejected) {
    const struct {
        std::string section, key, value;
    } cases[] = {
        {"molecules", "min_molecules_per_cell", "50.5"},
        {"segmentation", "n_clusters", "2147483648.0"},  // one past INT_MAX
        {"molecules", "min_molecules_per_cell", "1__000"},  // separators must sit between digits
    };
    for (const auto& c : cases) {
        EXPECT_THROW_MSG(load_config_text("[" + c.section + "]\n" + c.key + " = " + c.value + "\n"),
                         std::runtime_error,
                         "Invalid value '" + c.value + "' for config key '" + c.key + "': expected an integer");
    }
}
