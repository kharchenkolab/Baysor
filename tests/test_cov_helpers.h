// Shared helpers for the coverage test files (COV-*) and the BUG-* regression
// tests.
//
// The core of this header is portable: it uses no POSIX-only headers or
// calls (no getpid(), unistd.h, ...), so the test sources that include it
// stay buildable on Windows as well. Tests that are inherently POSIX
// (subprocess runs, setrlimit, /dev/full) guard themselves with
// #ifndef _WIN32; everything else relies on the helpers below. The CLI
// subprocess runner at the bottom is POSIX-only (sh-style quoting,
// WEXITSTATUS) and is compiled only for !defined(_WIN32) with BAYSOR_CLI_PATH
// defined; on Windows or without the CLI path the tests that use it skip.
#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
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

#include <spdlog/sinks/base_sink.h>
#include <spdlog/spdlog.h>

#include "baysor/utils/general.h"

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)
#include <sys/wait.h>
#endif

namespace baysor_test {

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Temporary directories
// ---------------------------------------------------------------------------

/// Create a fresh, uniquely named directory under the system temp directory.
/// Portable unique name: process-wide counter + timestamp + random suffix
/// (deliberately no getpid()).
inline fs::path make_unique_dir(const std::string& tag) {
    static std::atomic<int> counter{0};
    std::random_device rd;
    const auto ticks = std::chrono::steady_clock::now().time_since_epoch().count();
    const fs::path base = fs::temp_directory_path();
    for (int attempt = 0; attempt < 100; ++attempt) {
        fs::path dir = base / ("baysor_" + tag + "_" +
                               std::to_string(counter.fetch_add(1)) + "_" +
                               std::to_string(ticks) + "_" +
                               std::to_string(rd()));
        std::error_code ec;
        if (fs::create_directories(dir, ec)) return dir;
        if (ec && ec != std::errc::file_exists) continue;
        if (!ec && fs::is_directory(dir)) return dir;  // raced, but usable
    }
    throw std::runtime_error("test_cov_helpers: could not create a temp dir for tag '" +
                             tag + "'");
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
//
// Lives in a nested `cli` namespace on purpose: the COV-5 CLI tests keep
// their own local `run_cli`/`write_text` for TempDir arguments, and ADL on
// baysor_test::TempDir would make those calls ambiguous if these helpers
// sat directly in baysor_test.

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)
namespace cli {

/// Exit code plus captured stdout/stderr of a CLI subprocess run.
struct CliResult {
    int exit_code = -1;  // -1 = process did not exit normally (e.g. signal)
    std::string out;     // stdout
    std::string err;     // stderr
};

/// Read a whole file as binary text (empty string when unreadable).
inline std::string read_text_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

/// Write `content` to `tmp/name`; returns the file's full path.
inline std::string write_text(const TempDir& tmp, const std::string& name,
                              const std::string& content) {
    const fs::path p = tmp.path / name;
    std::ofstream f(p);
    f << content;
    return p.string();
}

/// Run the instrumented `baysor` binary (BAYSOR_CLI_PATH) with `args` in a
/// POSIX shell, capturing stdout/stderr into `tmp` and the exit code.
/// Shells out with sh-style quoting and decodes exit codes via WEXITSTATUS.
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
