# Revision history

Newest first.

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
