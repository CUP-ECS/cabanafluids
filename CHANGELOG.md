# Revision history

Newest first.

**2026-11-04 — tested the hypre solvers and fixed their stencil.** No test had
run the hypre pressure solve. `tstVelocityCorrector` now runs every case with
`Reference`/`Jacobi` and, when Cabana has hypre, with hypre's `PCG`/`Diagonal`,
`PCG`/`PFMG`, `GMRES`/`PFMG`, `BiCGSTAB`/`Diagonal`, `PFMG`, and `SMG`, and
checks that each hypre solve matches the reference one. Every hypre solver
failed: Cabana gives hypre the grid in (j, i) order but passes stencil offsets
through as (i, j), so the wall entries of the pressure matrix landed on the
wrong neighbors. CUP-ECS/Cabana PR #51 fixes this in the Cabana fork that
CabanaFluids builds against, and CI's image pins that fix. Hypre's `Jacobi` preconditioner stalls `PCG` and `GMRES` on the
all-solid test system, and `GMRES`/`Diagonal` reaches the `10 * sqrt(ncells)`
iteration cap, so the tests use the pairs above. When the example names a
hypre solver without `-c`, `createVelocityCorrector` now picks `PFMG` for the
Krylov solvers and none for `PFMG` and `SMG`; `-c` still wins, and the default
solver stays `Reference` with `Jacobi`. The example and the test main
now initialize hypre once after `MPI_Init` and finalize it before
`MPI_Finalize`; each hypre corrector used to initialize it, and nothing
finalized it.

**2026-11-04 — CI's image builds for a portable target.** Spack built each
image for the CPU of whichever runner built it: `zen2` for the old image and
`x86_64_v4`, which uses AVX-512, for the one PR #75 published. Test runners
without AVX-512 then crashed in the image's own `tar` during checkout, so no CI
job could build. `configs/github/spack.yaml` now requires `target=x86_64_v3`,
which every current GitHub-hosted x86 runner supports.

**2026-11-04 — CI's image pins the fixed Cabana.** The CUP-ECS Cabana fork
gave hypre its grid in reversed axis order but its stencil offsets unchanged,
so every hypre solver solved a system with i and j neighbors swapped; its PR
#51 fixes that. `configs/github/spack.yaml` now pins that merge commit, so the
CI image rebuilds with the fix rather than picking up `master` whenever it
next happens to rebuild.

**2026-11-04 — CI requires hypre.** CI's container builds Cabana with hypre,
but if it stopped doing so the hypre tests would compile out and CI would still
pass, and CI's log cannot show whether they ran. The new
`CabanaFluids_REQUIRE_HYPRE` option makes configure fail when Cabana lacks
hypre, and `configs/github/host-config.cmake` sets it. Local builds keep hypre
optional.

**2026-11-04 — Silo output records its time and makes its directories.**
`SiloWriter` passed the solver's `double` time where Silo's `DBOPT_TIME` reads
a `float`, so the `time` field held half of a double's bits, and it passed
`dt` as `DBOPT_DTIME`, so `dtime` held the step size. Both fields now hold the
time. The writer also wrote nothing unless `data/raw/` already existed and
still let the run exit 0; rank 0 now creates the directories, and a Silo file
that cannot be created throws `std::runtime_error` on every rank. The new
`SiloWriter` test suite covers both.

**2026-11-04 — updated the Spack environments and added a tutorial one.**
`configs/llnl-lassen/` is gone, because LLNL retired Lassen. The generic CUDA
environments in `configs/generic/` asked for upstream `cabana +hypre`, whose
CMake rejects hypre 3, and for the obsolete `hypre-cmake` package and
`concretization` key; they now take the CUP-ECS Cabana fork from the `cupecs`
repo, plain `hypre`, Kokkos 4 with C++17, and Silo 4.12, as
`configs/github/spack.yaml` does. `configs/tutorial/spack.yaml` replaces the
2022 UAO tutorial environment on `main` with a CPU-only one that students can
build on Linux or macOS, using Kokkos' Serial and Threads backends. Many
tutorial machines have GPUs, so `configs/tutorial-cuda/spack.yaml` is its CUDA
counterpart, for Volta by default.

**2026-11-04 — the example runs on Kokkos' Threads backend.** The unit tests
already built for Threads, but `createSolver` rejected `-p threads` as an
invalid backend, so the example could not run on it. It now builds a
`Kokkos::Threads` solver under the same guards as `serial` and `openmp`, and
`tstSolver` checks that it does. The `-p` help text and `README.md` list the
backend names, and `README.md` now gives the real default, `default`, which
picks Kokkos' default execution space and so runs Threads in a build that
enables it.

**2026-11-04 — `actions/checkout@v7` in every workflow.** The build-and-test
and `clang-format` workflows used `actions/checkout@v3`, which ran on a Node
release GitHub has retired, while the container workflow already used v7. All
three now use v7.

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
