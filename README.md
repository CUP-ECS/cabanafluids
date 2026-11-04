# CabanaFluids
This directory contains code for a relatively simple finite difference 
fluid advection solver for exploring communication issues on modern architectures
(particularly GPUs). The main goal is to look at different neighbor collective
and GPU communication approaches.

Computationally, the benchmark advects a material feature (that doesn't otherwise affect
fluid flow, e.g. by changing pressure) using the incompressible Euler fluid flow equations. 
Consider, for example, something like a dye being carried through a tank of water or a 
fragrance wafting across a room. 

The main elements of the benchmark are:
  * Solution of the pressure gradient at each timestep to maintain 
    incompressibility on a staggered marker-and-cell (MAC) grid. The
    benchmark has two implementations: (1) the HYPRE structured solvers
    through `Cabana::Grid`, when Cabana is built with HYPRE, or (2) the
    `Cabana::Grid` reference preconditioned CG solver (the default), in
    which different MPI approaches for the halo exchange are explored.
  * Semi-Lagrangian advection: 3rd-order Runge-Kutta back-traces through
    linearly interpolated velocities, and cubic B-spline interpolation of
    the advected quantity.
  * A time step limited each step so that no back-trace travels more than
    one cell; the command-line time step is a maximum.
  * Free-slip walls, with ghost cells outside the domain filled by
    reflection after every halo exchange.

Sources:
  - Fluid Simulation for Computer Graphics by Bridson
  - Incremental Fluids in Kokkos (git@github.com:pkestene/incremental-fluids-kokkos.git)

## Building and Testing

CabanaFluids requires CMake 3.18 or later, a C++17 compiler, MPI, Cabana
built with MPI and grid support, and Silo. HYPRE is optional; configure
with `-DCabanaFluids_REQUIRE_HYPRE=ON` to fail when Cabana lacks it. The unit
tests need GoogleTest; configure with `-DCabanaFluids_ENABLE_TESTING=OFF`
to skip them.
The `configs` directory holds Spack environments for several systems.

```
cmake -B build
cmake --build build
ctest --test-dir build
```

CTest runs each test suite at 1, 2, and 4 MPI ranks on every enabled Kokkos
backend. Set
`MPIEXEC_PREFLAGS` if your MPI needs extra arguments, as
`configs/github/host-config.cmake` does for OpenMPI in CI.

## Running

`build/examples/advection` runs a 2D problem in a unit box with an inflow
source that injects material moving in the +x direction. Run it with `--help` for the full option list. Common
options:

| Option | Meaning | Default |
| --- | --- | --- |
| `-p` | Kokkos backend: `serial`, `openmp`, `cuda`, or `hip` | `serial` (`cuda` in CUDA builds) |
| `-n` | Cells per dimension | 128 |
| `-t` | Simulated time | 4.0 |
| `-i` | Maximum time step | 0.005 |
| `-f` | Output frequency in steps | 20 |
| `-m`, `-c` | Pressure solver and preconditioner | `Reference`, `Jacobi` |

The solver writes Silo files every `-f` steps for viewing in VisIt.

## Status and Future Work

The advected quantity's cubic B-spline interpolation does not reproduce
grid values, so every step adds numerical diffusion. `FUTURE.md` records
this bug, the bugs already fixed, a plan for the remaining unit and
convergence tests, and a plan to replace semi-Lagrangian advection with a
conservative MUSCL scheme and SSP-RK3 time integration.
