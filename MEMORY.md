# Session memory

Reasoning worth carrying between agent sessions that belongs in no other
file. Written 2026-11-01, updated 2026-11-04. Verify before acting on
it; the files it points at win where the two disagree.

## Current state

- 2026-11-04: work toward the first release follows Section 0 of
  `FUTURE.md`, which records each step's status and pull request. Steps
  run in parallel on separate branches from `develop`, so their
  `CHANGELOG.md` entries conflict at the top as they merge. Issue #54
  and the numerics wait until after the release.

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
  upstream Cabana's CMake rejects hypre 3. Since 2026-11-04 the local
  Spack environment builds the same fork with `+hypre +threads`, so
  local builds compile the hypre path too.
- The tests run with Kokkos' Serial backend locally and in CI, and with
  Threads locally. OpenMP, CUDA, and HIP builds of the harness are
  untested.

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
  `../.spack-env/view`, with its `mpicxx` and `mpicc` as the compilers
  and `-DCMAKE_PREFIX_PATH` set to the view; that suffices since the
  2026-11-04 rebuild. Use a build directory outside the repository. The
  older recipe copied `CMAKE_PREFIX_PATH` from the Spack build's
  `CMakeCache.txt` under `build-darwin-*/`. Sourcing that build's `spack-build-env.txt`
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
  `develop`. It suffices for the 64x64 test mesh with every solver PR
  #70 tests except `GMRES` with `Diagonal`; nobody has checked it at
  production sizes.
- `FUTURE.md` holds the remaining test plan and the advection rewrite.
- Waiting on the author: the release version number.
- Decided 2026-11-04: the example keeps Cabana's built-in CG
  (`Reference`) as its default solver; the second local backend is
  Threads, because Apple Clang lacks OpenMP; Lassen and Quartz are
  retired, so their configurations go; a tutorial Spack environment
  replaces the UAO one from `main`.
- Decided 2026-11-04, later: Cabana bugs are fixed in the CUP-ECS fork,
  which every configuration uses, and the author upstreams them; the
  hypre default preconditioner changes; CI's configure requires hypre;
  tutorials get CPU and CUDA environments, Volta by default, pinned to
  `develop` until the release; MPI tests need not run Threads with
  more than one thread; roundoff-level differences between backends are
  acceptable; the Silo output bugs are fixed before the release.
