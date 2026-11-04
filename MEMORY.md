# Session memory

Reasoning worth carrying between agent sessions that belongs in no other
file. Written 2026-11-01, updated 2026-11-04. Verify before acting on
it; the files it points at win where the two disagree.

## Current state

- 2026-11-04: work toward the first release started on
  `prepare-release`, branched from `origin/develop`. The plan is
  Section 0 of `FUTURE.md`: reconcile `main`, test hypre locally and in
  CI, remove the CajitaFluids name, and prepare the release. Nothing
  beyond the plan is done yet. Issue #54 and the numerics wait until
  after the release.

- `cabanafluids-cleanup` was rebased onto `develop` and merged as PR
  #60; its branches are deleted. `CHANGELOG.md` records what changed
  and why.
- CI runs in `ghcr.io/cup-ecs/cabanafluids-tester:latest`, rebuilt from
  `configs/github/spack.yaml` by PRs #61 and #62 after the old
  `cajitafluids-tester` image broke `actions/checkout`. Pull requests
  only build the image; a push to `develop` publishes it, which takes
  about an hour. Changing the workflow or the image needs a human.
- The image builds the CUP-ECS Cabana fork's `master` with `+hypre`
  from the `cupecs` repo in COMPASS-PSAAP/spack-packages, because
  upstream Cabana's CMake rejects hypre 3. Local builds use the same
  fork without hypre, so only CI compiles the hypre solver path, and no
  test runs it: the `VelocityCorrector` tests ask for `"Reference"`.
- The tests have run only with Kokkos' Serial backend, on macOS and in
  CI's Linux container, where all 18 pass. OpenMP, CUDA, and HIP builds
  of the harness are untested.

## Design decisions

- The rebase kept `develop`'s CMake build and Cabana test harness over
  this branch's own BLT replacement, because `develop` owns the testing
  infrastructure.
- Tests keep their typed GoogleTest fixtures, with `MeshDeviceTypes` in
  `tests/tstDriver.hpp` reduced to the one device the harness compiles
  for. Converting to Cabana's `TEST( TEST_CATEGORY, ... )` style would
  rewrite every test for no behavioral gain.
- The vendored `cmake/test_harness/test_harness.cmake` is edited to run
  exactly 1, 2, and 4 ranks. Do not restore Cabana's
  `MPIEXEC_MAX_NUMPROCS` logic; it drops 4-rank runs on small machines
  and adds runs at odd rank counts the tests do not expect.

## Procedures easy to get wrong

- `rerere` is enabled in this clone. A trial rebase resolved by "take one
  side" records those resolutions and replays them silently later; clear
  the `.git/rr-cache` entries a trial created.
- Local build: configure against the Spack view in
  `../.spack-env/view`, with `mpicxx` as the compiler and
  `CMAKE_PREFIX_PATH` copied from the Spack build's `CMakeCache.txt`
  under `build-darwin-*/`. Sourcing that build's `spack-build-env.txt`
  aborts the shell, and the file holds an API key; never commit or quote
  it.
- A new test suite is `tests/tst<Name>.hpp` plus an entry in `MPI_TESTS`
  in `tests/CMakeLists.txt`; shared fixtures live in `*Fixture.hpp`,
  because the harness claims every `tst<Name>.hpp` name.
- `ctest -N` reporting zero tests means the harness found no devices.
  `CABANA_SUPPORTED_DEVICES` is defined in the root `CMakeLists.txt`
  because Cabana does not export it.
- Local `clang-format` is version 23; CI uses 14.
- Check a `spack.yaml` change with fresh repo clones, user config off
  (`SPACK_USER_CONFIG_PATH`, `SPACK_DISABLE_LOCAL_CONFIG=1`), and
  `spack concretize --fresh`. Spack checks out a pinned repo `commit`
  only on a fresh clone, and local installs mask the versions CI gets.
  A clean concretize still does not prove the container builds.

## Open questions

- `VelocityCorrector` caps CG at `10 * sqrt(ncells)` iterations, from
  `develop`. It suffices for the 64x64 test mesh; nobody has checked it
  at production sizes or with HYPRE solvers.
- A `VelocityCorrector` test that selects a hypre solver would exercise
  the hypre path in CI; it needs `Cabana_ENABLE_HYPRE` guards, since
  local builds lack hypre.
- `FUTURE.md` holds the remaining test plan and the advection rewrite.
- Waiting on the author: the release version number.
- Decided 2026-11-04: the example keeps Cabana's built-in CG
  (`Reference`) as its default solver; the second local backend is
  Threads, because Apple Clang lacks OpenMP; Lassen and Quartz are
  retired, so their configurations go; a tutorial Spack environment
  replaces the UAO one from `main`.
