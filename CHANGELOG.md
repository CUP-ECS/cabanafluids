# Revision history

Newest first.

**2026-11-04 — updated the Spack environments and added a tutorial one.**
`configs/llnl-lassen/` is gone, because LLNL retired Lassen. The generic CUDA
environments in `configs/generic/` asked for upstream `cabana +hypre`, whose
CMake rejects hypre 3, and for the obsolete `hypre-cmake` package and
`concretization` key; they now take the CUP-ECS Cabana fork from the `cupecs`
repo, plain `hypre`, Kokkos 4 with C++17, and Silo 4.12, as
`configs/github/spack.yaml` does. `configs/tutorial/spack.yaml` replaces the
2022 UAO tutorial environment on `main` with a CPU-only one that students can
build on Linux or macOS, using Kokkos' Serial and Threads backends.

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
