# Revision history

Newest first.

**2026-11-04 — merged `main` into `develop`.** `main` held one commit
`develop` lacked, `d5a0b83` from 2022, so the release merge of `develop` into
`main` conflicted in `ProblemManager.hpp`, `SiloWriter.hpp`, and `Solver.hpp`.
The merge keeps `develop`'s `Cabana::Grid` code and ports `main`'s fixes: the
I- and J-face velocity getters return the face-array view types, and Silo
output files end in `.silo` instead of `.pdb`, which changes the example's
output file names. `develop` already balanced the profiling regions in
`Solver`, so its version stays. The merge drops
`configs/tutorial-uao-cuda/spack.yaml`; the `update-spack-configs` branch
replaces it with `configs/tutorial/`.

**2026-11-01 — rebased onto `develop` and its test harness.** `develop`
replaced BLT with plain CMake and Cabana's test harness, as this branch had
done separately; the rebase keeps `develop`'s build. The harness registered no
tests, because Cabana does not export the device list it loops over, and
CTest passed with nothing to run. `CMakeLists.txt` now defines that list, each
suite moved from `tst<Name>.cpp` to the `tst<Name>.hpp` the harness compiles
once per Kokkos backend, and the shared fixtures moved to `*Fixture.hpp`.
Every suite runs at 1, 2, and 4 ranks regardless of core count, so the
4-rank runs that exposed earlier decomposition bugs never drop out. The
pressure solve's iteration cap from `develop` counted cells in one dimension
only, which stopped CG before it converged on the test mesh; it now counts
all cells.

**2026-11-01 — inflow as an overlap-weighted rate.** `InflowSource` tested
only cell and face centers against its box, so a small box could miss the mesh
entirely, and it forced a minimum quantity rather than adding a rate. Each cell
and face now receives a share in proportion to how much of its control volume
the box covers, so every box of positive area reaches at least one cell and
the total injected is `rate * dt * area` at any resolution. The example's `-q`
option is now a rate per unit area per unit time, which changes its output.
Uniform gravity turned out to have no effect on a constant-density fluid in a
closed box; `FUTURE.md` records this, and `tstVelocityCorrector` tests it.

**2026-10-30 — adopted `COMPASS-PSAAP/compass-template`.** Added `AGENTS.md` with
`CLAUDE.md` as a symlink, this changelog, the BSD 3-Clause `LICENSE`, the
pull-request and issue templates, and placeholder Spack environments under
`spack-envs/llnl/`. The template's `clang-format` workflow replaces
`clang-format-check.yml`, so CI now checks the whole repository with
`clang-format` 14 instead of `src`, `tests`, and `examples` with version 13.
`.gitignore` follows the template and keeps `build-*` for Spack build symlinks.
This repository keeps its own LLVM-based Allman `.clang-format`, because the
template's style would reformat every source file. The existing `configs/`
directory, the build-and-test workflow, and the CI-container workflow stay.
