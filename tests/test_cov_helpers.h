// Shared test helpers. Everything except the CLI subprocess runner at the
// bottom (POSIX only, needs BAYSOR_CLI_PATH) is portable.
#pragma once

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <spdlog/sinks/base_sink.h>
#include <spdlog/spdlog.h>

#include "baysor/utils/general.h"
#include "baysor/utils/thread_pool.h"

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)
#include <sys/wait.h>
#endif

/// Expects `statement` to throw `exception` whose what() contains `substr`.
#define EXPECT_THROW_MSG(statement, exception, substr)                            \
    EXPECT_THROW(                                                                 \
        try { statement; } catch (const exception& e_) {                          \
            EXPECT_NE(std::string(e_.what()).find(substr), std::string::npos)     \
                << "message: " << e_.what();                                      \
            throw;                                                                \
        },                                                                        \
        exception)

namespace baysor_test {

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Files
// ---------------------------------------------------------------------------

/// Create a fresh, uniquely named directory under the system temp directory.
inline fs::path make_unique_dir(const std::string& tag) {
    std::random_device rd;
    for (;;) {
        const fs::path dir = fs::temp_directory_path() /
            ("baysor_" + tag + "_" + std::to_string(rd()) + std::to_string(rd()));
        if (fs::create_directories(dir)) return dir;
    }
}

/// Whole file as a string ("" when unreadable).
inline std::string read_text_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

/// Values of column `name` of a simple (unquoted) CSV file.
inline std::vector<std::string> csv_column(const fs::path& p, const std::string& name) {
    auto split = [](const std::string& line) {
        std::vector<std::string> fields;
        std::stringstream ss(line);
        for (std::string field; std::getline(ss, field, ',');) fields.push_back(field);
        return fields;
    };
    std::ifstream f(p);
    std::string line;
    std::getline(f, line);
    const auto header = split(line);
    const size_t col = std::find(header.begin(), header.end(), name) - header.begin();
    if (col == header.size()) throw std::runtime_error("no column '" + name + "' in " + p.string());
    std::vector<std::string> values;
    while (std::getline(f, line)) {
        const auto fields = split(line);
        values.push_back(col < fields.size() ? fields[col] : "");
    }
    return values;
}

/// Per-test unique directory, removed on destruction (RAII).
class TempDir {
public:
    fs::path path;

    explicit TempDir(const std::string& tag = "test") : path(make_unique_dir(tag)) {}
    ~TempDir() {
        std::error_code ec;
        fs::remove_all(path, ec);
    }

    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;

    /// Path of a file inside the directory, as a string.
    std::string file(const std::string& name) const { return (path / name).string(); }

    /// Write `content` to the file `name` inside the directory; returns its path.
    std::string write(const std::string& name, const std::string& content) const {
        std::ofstream(path / name) << content;
        return file(name);
    }
};

/// Sets the thread-pool size for the guard's lifetime.
class PoolSizeGuard {
public:
    explicit PoolSizeGuard(int n) : old_(baysor::thread_pool_size()) { baysor::set_thread_pool_size(n); }
    ~PoolSizeGuard() { baysor::set_thread_pool_size(old_); }
    PoolSizeGuard(const PoolSizeGuard&) = delete;
    PoolSizeGuard& operator=(const PoolSizeGuard&) = delete;

private:
    int old_;
};

// ---------------------------------------------------------------------------
// Global RNG isolation
// ---------------------------------------------------------------------------

/// Restores the project's default global RNG stream (reset_global_xoshiro_rng
/// with its default seed) on construction *and* destruction. Tests that
/// reseed or advance the global xoshiro hold one of these so their results do
/// not depend on which tests ran before them.
struct GlobalRngGuard {
    GlobalRngGuard() { baysor::reset_global_xoshiro_rng(); }
    ~GlobalRngGuard() { baysor::reset_global_xoshiro_rng(); }

    GlobalRngGuard(const GlobalRngGuard&) = delete;
    GlobalRngGuard& operator=(const GlobalRngGuard&) = delete;
};

// ---------------------------------------------------------------------------
// Log capture
// ---------------------------------------------------------------------------

/// spdlog sink that records every payload (level and metadata stripped) so
/// tests can assert on the exact text of warnings/infos.
class CapturingSink : public spdlog::sinks::base_sink<std::mutex> {
protected:
    void sink_it_(const spdlog::details::log_msg& msg) override {
        data_.append(msg.payload.data(), msg.payload.size());
        data_.push_back('\n');
    }
    void flush_() override {}

public:
    std::string data() const { return data_; }
    void clear() { data_.clear(); }

private:
    std::string data_;
};

/// Swaps the spdlog default logger for one that feeds `sink`, and restores
/// the original logger on destruction (even when expectations fail).
struct LoggerGuard {
    std::shared_ptr<spdlog::logger> original;

    explicit LoggerGuard(std::shared_ptr<spdlog::sinks::sink> sink)
        : original(spdlog::default_logger()) {
        sink->set_level(spdlog::level::trace);
        auto logger = std::make_shared<spdlog::logger>("cov-capture", std::move(sink));
        logger->set_level(spdlog::level::trace);
        spdlog::set_default_logger(logger);
    }
    ~LoggerGuard() { spdlog::set_default_logger(original); }

    LoggerGuard(const LoggerGuard&) = delete;
    LoggerGuard& operator=(const LoggerGuard&) = delete;
};

// ---------------------------------------------------------------------------
// CLI subprocess runner (POSIX + BAYSOR_CLI_PATH only)
// ---------------------------------------------------------------------------

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)
namespace cli {

using baysor_test::read_text_file;

/// Exit code plus captured stdout/stderr of a CLI subprocess run.
struct CliResult {
    int exit_code = -1;  // -1 = process did not exit normally (e.g. signal)
    std::string out;     // stdout
    std::string err;     // stderr
};

inline std::string write_text(const TempDir& tmp, const std::string& name,
                              const std::string& content) {
    return tmp.write(name, content);
}

/// Run the `baysor` binary (BAYSOR_CLI_PATH) with `args` in a POSIX shell,
/// capturing stdout/stderr into `tmp` and the exit code.
inline CliResult run_cli(const TempDir& tmp, const std::string& args) {
    const fs::path out_p = tmp.path / "stdout.txt";
    const fs::path err_p = tmp.path / "stderr.txt";
    const std::string cmd = "'" + std::string(BAYSOR_CLI_PATH) + "' " + args +
                            " > '" + out_p.string() + "' 2> '" +
                            err_p.string() + "'";
    const int status = std::system(cmd.c_str());

    CliResult r;
    if (status >= 0 && WIFEXITED(status)) {
        r.exit_code = WEXITSTATUS(status);
    }
    r.out = read_text_file(out_p);
    r.err = read_text_file(err_p);
    return r;
}

}  // namespace cli
#endif  // !defined(_WIN32) && defined(BAYSOR_CLI_PATH)

}  // namespace baysor_test
