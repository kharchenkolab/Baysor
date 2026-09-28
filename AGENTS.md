# Baysor workspace

Paths below are relative to the active checkout or worktree.

- Work within this Baysor checkout or its assigned worktree.
- Read `README.md` and relevant `docs/` before choosing build or test commands.
- Build configuration lives in `CMakeLists.txt`, `CMakePresets.json` and
  `configure.sh`; inspect the relevant configuration before running checks.
- Source code is in `src/` and `include/`; tests are in `tests/`.
- Preserve existing local build directories and unrelated changes.
- Do not assume a dispatcher or PLAN.md exists. Apply a dispatched-job protocol
  only when a dispatcher explicitly supplies a job.
