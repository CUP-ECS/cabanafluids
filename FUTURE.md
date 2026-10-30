# Future Work for CabanaFluids

This document records known bugs, a plan for unit tests, and a plan for
moving the advection step to a robust, well-documented finite difference
scheme. Line numbers refer to commit `b5284b9`.

## 1. Known and Suspected Bugs

Fix each bug red-green: write the failing test from Section 2 first, then
fix the code.

### Correctness bugs in the solver

1. **J-face boundary condition applied to `u`.**
   `src/VelocityCorrector.hpp:264` calls `bc( FaceJ(), u, ... )` instead of
   `bc( FaceJ(), v, ... )`. After the pressure correction the code never
   zeroes `v` on the top and bottom walls, so fluid leaks through them.
   It also zeroes entries of `u` at J-face indices, including a ghost row.

2. **RK3 back-trace is only second order.**
   `src/TimeIntegrator.hpp:58-59` computes the third stage position from
   `v0`. Ralston's RK3, whose weights (2/9, 3/9, 4/9) the code uses,
   needs `x2 = x0 - 0.75 * dt * v1`. With `v0` the order condition
   `b3 * a32 * c2 = 1/6` fails. Line 61 has a second error in the 3D path:
   it writes `x1[2]` with factor 0.5 instead of `x2[2]` with 0.75.

3. **Cubic interpolation smooths the field every step.**
   `src/TimeIntegrator.hpp:114` interpolates the advected quantity with
   `Cabana::Grid` order-3 splines. These are cubic B-splines, which
   approximate rather than interpolate: they do not reproduce grid values.
   Even with zero velocity, each step blurs the field with weights
   (1/6, 4/6, 1/6) per dimension. This adds large numerical diffusion.

4. **Ghost cells at physical walls are never filled.**
   Halo gathers fill ghosts only between ranks. Ghosts outside the
   physical domain keep their initial value of zero. Back-traces near a
   wall interpolate those zeros, so the quantity and the velocity decay
   along walls.

5. **Back-traces can leave the allocated halo.**
   `src/TimeIntegrator.hpp:100-116` never clamps the traced point.
   `src/Solver.hpp:95-105` picks `dt` once, from the inflow velocity and
   body force, but the pressure correction can produce faster flow. If a
   trace moves more than about one cell, the cubic stencil reads past the
   3-cell halo: an out-of-bounds view access, which corrupts memory or
   crashes on a GPU. The time step should be recomputed from the actual
   maximum velocity (a global `MPI_Allreduce`) each step.

6. **Singular pressure system with all-solid walls.**
   The example sets all four walls `SOLID`
   (`examples/advection.cpp:447-449`), so the matrix is a pure Neumann
   Laplacian with a constant null space. CG converges only if the
   right-hand side has zero mean. Round-off breaks that, and nothing pins
   the pressure or removes the mean.

### Bugs in infrastructure and the driver

7. **Unbalanced profiling regions.** In `src/TimeIntegrator.hpp:164-177`,
   the `popRegion` for "Advance" sits inside the 3D-only branch. In 2D the
   final `popRegion` closes "Advance" and leaves "TimeIntegrator::Step"
   open every step.
8. **HIP backend does not compile.** `src/Solver.hpp:337` has
   `Kokkos : Experimental::HIP` (one colon).
9. **3D advection call is missing an argument.**
   `src/TimeIntegrator.hpp:158` omits `bc`. This is latent until 3D is
   built.
10. **Mesh divisibility check uses an absolute tolerance.**
    `src/Mesh.hpp:61-63` compares extents against `10 * epsilon`, which
    does not scale with the domain size. Large domains can throw a
    spurious `logic_error`. Use a relative tolerance.
11. **Initial output is overwritten.** `src/Solver.hpp:155` writes step 0,
    then the loop at lines 171-174 writes step 0 again after the first
    step. The `strdup( "Mesh" )` calls leak.
12. **Velocity output reads stale ghosts.** `_applyPressure` updates owned
    faces but does not gather afterward. `src/SiloWriter.hpp:183-194` then
    interpolates cell-centered velocity from ghost faces on rank
    boundaries that predate the correction.
13. **Initializer receives local indices.** `src/ProblemManager.hpp:214`
    passes local `(i, j)` to the init functor as `coords`. A functor that
    expects global indices behaves differently on each rank. Pass global
    indices or document the contract.
14. **Inconsistent index conversion.**
    `src/VelocityCorrector.hpp:228-229` builds a `Cell` L2G and applies it
    to face indices; `src/Solver.hpp:221` uses face L2G. The two agree for
    non-periodic grids only.
15. **Command-line parsing in `examples/advection.cpp`.**
    - Line 214 validates `global_bounding_box[3]` before assigning it;
      it should check `[2]`.
    - Line 241 validates `t_final` instead of `delta_t`.
    - `-w` sets the inflow width, so the write frequency cannot be set;
      `-h` sets inflow height, so only `--help` prints help.
    - The help text lists `-t` twice, and several stated defaults differ
      from the code (inflow width and velocity).
    - The parser accepts `pthreads`, which `createSolver` rejects, and
      rejects `hip`, which `createSolver` supports.
16. **Documentation.** The README title says "Poisson MPI Benchmark", and
    the example header says the domain is open in X, but the driver closes
    all four walls.

### Bugs in the existing tests

17. `tests/tstMesh.cpp` `LocalGridSetup` uses `dimNumBlock( 1 )` for
    dimension 0 and assumes the cell count divides evenly among ranks.
18. `tests/tstProblemManager.cpp` `HaloTest` encodes each rank's *local*
    indices in the values and compares them against the receiver's local
    indices. The local indices of a shared cell differ between ranks, so
    the test passes only on one rank, where there are no neighbors.
19. All four tests in `tests/tstBoundaryConditions.cpp` are empty, and
    CTest runs every test on a single rank.

## 2. Unit Test Plan

### Infrastructure

- Register each test at 1, 2, and 4 ranks, so halo and decomposition
  code runs under real communication. The build uses plain CMake and
  CTest, so extend the loop in `tests/CMakeLists.txt` with an inner loop
  over rank counts that calls
  `add_test( NAME ${_test}Tests_${n} COMMAND ${MPIEXEC_EXECUTABLE}
  ${MPIEXEC_NUMPROC_FLAG} ${n} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:tst${_test}>
  ${MPIEXEC_POSTFLAGS} )` and sets the test's `PROCESSORS` property to
  `n`. `configs/github/host-config.cmake` already sets `MPIEXEC_PREFLAGS`
  so OpenMPI runs as root and oversubscribes in the CI container, and the
  GitHub workflow runs every registered test through `ctest`.
- Add a small helper header for tests: copy a `Cabana::Grid` array to the
  host, compute a global max, sum, or L2 norm over owned entities with
  `MPI_Allreduce`, and fill an array from a function of global position.
- Encode global, not local, indices in any value that crosses ranks.
- Keep test meshes small (32x32 to 128x128) so the suite finishes in
  seconds on CPU backends.

### Tests by component

**Mesh**
- Cell size, global index bounds, and owned extents for rank counts that
  give uneven decompositions.
- A large domain (for example `[0, 1e6]`) constructs without throwing;
  a non-divisible domain throws (bug 10).

**ProblemManager**
- `advance` swaps current and next for every field, not only the cell
  field.
- Halo exchange: write `f(global i, global j)` into owned entities,
  gather, and check every ghost entity that has a neighbor (bug 18).
- The init functor receives the cell-center and face-center coordinates
  that `local_mesh.coordinates` reports (bug 13).

**BoundaryCondition**
- Velocity: with `SOLID` walls, every wall-normal face velocity is zero
  after the functor runs; with `FREE` walls, values are unchanged. Run
  separate cases for I-faces and J-faces (bug 1).
- Matrix: check the diagonal and off-diagonal entries for interior, edge,
  and corner cells under each wall type. A solid edge cell has diagonal
  `3 * scale`, a solid corner `2 * scale`, and a free edge keeps
  `4 * scale` with the outside neighbor removed.

**Interpolation**
- Linear interpolation reproduces any linear field exactly.
- The chosen higher-order interpolant reproduces grid values at grid
  points; a failing version of this test documents bug 3.
- Face-centered interpolation lands on the right staggered locations.

**TimeIntegrator**
- `rk3` with constant velocity returns `x0 - dt * u` exactly.
- `rk3` in the solid-body rotation `u = (-y, x)`: the trace error falls by
  about 8x when `dt` halves (bug 2).
- `advect` with zero velocity leaves the field unchanged to round-off
  (bug 3).
- `advect` near a wall, and with a velocity that would trace outside the
  halo, stays in bounds and produces finite values (bugs 4 and 5). Run
  this under a Kokkos debug build with bounds checking.

**VelocityCorrector**
- After `correctVelocity`, the discrete divergence of every owned cell is
  below the solver tolerance, and every solid-wall normal velocity is zero.
- Projection is idempotent: a second call changes nothing beyond the
  solver tolerance.
- Given `u = u_df + grad(phi)` for a known divergence-free `u_df`,
  projection recovers `u_df` to second order in `h`.
- With all-solid walls and a right-hand side that has a small nonzero mean,
  the solve still converges (bug 6).

**Solver (integration tests)**
- Closed box, no inflow, no body force: total quantity is conserved to the
  scheme's tolerance, and the quantity stays within its initial min and
  max.
- One-rank and four-rank runs of the same short problem agree to round-off.
- The time step adapts so the CFL number never exceeds the configured
  limit.
- Command-line parsing: move `parseInput` into a testable function and
  check the cases in bug 15.

### Convergence tests

Give these tests the CTest label `convergence` (the `LABELS` test
property) and exclude them from CI with `ctest -LE convergence`.

- Translate a smooth Gaussian with constant velocity and compare to the
  exact solution. Measure the L1 error at 32, 64, and 128 cells and check
  the observed order.
- Rotate a slotted disk (Zalesak's test) for one revolution and check that
  the solution stays within bounds and the shape survives.

## 3. Plan for a Robust Finite Difference Scheme

### Recommendation

Keep the staggered marker-and-cell (MAC) grid and the Chorin pressure
projection; both are standard and already in place. Replace
semi-Lagrangian advection with a conservative flux-form finite difference
scheme:

- **Spatial discretization:** second-order upwind reconstruction (MUSCL)
  with a slope limiter. Use the van Leer or monotonized-central limiter.
  The limiter keeps the scheme total-variation diminishing, so it creates
  no new extrema, and it falls back to first-order upwind near steep
  gradients.
- **Time integration:** the three-stage strong-stability-preserving
  Runge-Kutta method of Shu and Osher (SSP-RK3). Each stage is a forward
  Euler step, and the stages combine convexly, so the limiter's
  no-new-extrema property survives the full step.
- **Time step:** choose `dt = C * h / max|u|` with `C <= 0.8`, computing
  `max|u|` each step with a global reduction.

### Why this scheme

- Conservation: flux form moves material between cells without creating
  or destroying it, which the integration tests can check exactly.
- Bounded solutions: the limiter and SSP time stepping keep quantities
  within their initial range, without the B-spline smoothing of the
  current code.
- Fixed stencils: each update reads two cells in each direction. Two halo
  cells always suffice, so traces cannot run off the allocated halo.
- Teachability: every piece fits in a few lines of code and appears in
  standard textbooks. The code can cite the equation each line implements.
- Benchmark value: SSP-RK3 needs one halo exchange per stage, three per
  step, with a fixed, predictable message pattern. This suits the
  project's goal of studying halo-exchange strategies.

The main cost is the CFL limit: semi-Lagrangian advection tolerates large
time steps, and the flux form does not. For this benchmark, predictable
stencils and correctness matter more than large steps.

### Discretization summary

For a cell-centered scalar `q` and face velocities `u` and `v`, one forward
Euler stage is

```
q_new(i,j) = q(i,j) - dt/h * ( F(i+1/2,j) - F(i-1/2,j)
                             + G(i,j+1/2) - G(i,j-1/2) )
```

The flux through the I-face at `i-1/2` takes the upwind reconstruction:

```
F(i-1/2,j) = u(i,j) * ( u(i,j) > 0 ? qL : qR )
qL = q(i-1,j) + 0.5 * phi(r_L) * ( q(i,j)   - q(i-1,j) )
qR = q(i,j)   - 0.5 * phi(r_R) * ( q(i+1,j) - q(i,j)   )
```

Here `u(i,j)` is the stored I-face velocity at `i-1/2`, `phi` is the
limiter, and `r` is the ratio of consecutive differences. The velocity
components use the same form on their own staggered control volumes, with
transport velocities averaged to the faces of those volumes.

SSP-RK3 combines three such stages:

```
q1 = q  + dt * L(q)
q2 = 3/4 q + 1/4 ( q1 + dt * L(q1) )
q_new = 1/3 q + 2/3 ( q2 + dt * L(q2) )
```

### Implementation steps

1. **Fix the bugs in Section 1** and land the unit tests, so the new scheme
   starts from a tested base.
2. **Fill physical ghost cells.** Add a `BoundaryCondition::fillGhosts`
   step after each halo gather: reflect normal velocity at solid walls,
   copy tangential velocity (free slip) or negate it (no slip), and apply
   zero-gradient extrapolation for the scalar at solid and free walls.
   Test it independently.
3. **Add an adaptive time step.** Compute `max|u|` with a Kokkos reduction
   and `MPI_Allreduce` each step. Remove the one-time estimate in
   `Solver`.
4. **Add flux-form advection** in a new `src/Advection.hpp`: a limiter
   function, a face-flux kernel for each direction, and a flux-divergence
   kernel. Document each function with the equation it implements and a
   textbook reference.
5. **Replace `TimeIntegrator::step`** with SSP-RK3. Store one extra stage
   array per field in `ProblemManager`, and gather and fill ghosts before
   each stage.
6. **Handle the pressure null space.** Subtract the mean of the right-hand
   side before the solve when every wall is solid, or pin one pressure
   value.
7. **Reduce the halo to two cells** once nothing needs three, and update
   the comment in `ProblemManager.hpp` that explains the halo depth.
8. **Keep the old scheme for comparison** only if benchmark users need it.
   Otherwise remove it, along with the spline interpolation of advected
   quantities. Keep linear interpolation for output.
9. **Write `docs/numerics.md`:** the governing equations, the MAC grid
   layout with index conventions, the scheme above, boundary handling, and
   the CFL condition, with each section pointing to the code that
   implements it.

### References

- F. H. Harlow and J. E. Welch, "Numerical calculation of time-dependent
  viscous incompressible flow of fluid with free surface," *Physics of
  Fluids* 8 (1965). The MAC grid.
- A. J. Chorin, "Numerical solution of the Navier-Stokes equations,"
  *Mathematics of Computation* 22 (1968). Pressure projection.
- B. van Leer, "Towards the ultimate conservative difference scheme V,"
  *Journal of Computational Physics* 32 (1979). MUSCL reconstruction.
- C.-W. Shu and S. Osher, "Efficient implementation of essentially
  non-oscillatory shock-capturing schemes," *Journal of Computational
  Physics* 77 (1988). SSP-RK3.
- R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*,
  Cambridge University Press (2002). Limiters and TVD theory, chapters 6
  and 19.
- R. Bridson, *Fluid Simulation for Computer Graphics*, 2nd ed., CRC Press
  (2015). The MAC grid and projection as used in this code.
