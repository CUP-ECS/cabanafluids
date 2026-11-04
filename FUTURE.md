# Future Work for CabanaFluids

This document records known bugs, a plan for unit tests, and a plan for
moving the advection step to a robust, well-documented finite difference
scheme.

## 1. Known Bugs

**Cubic interpolation smooths the field every step.**
`src/TimeIntegrator.hpp:155` interpolates the advected quantity with
`Cabana::Grid` order-3 splines. These are cubic B-splines, which
approximate rather than interpolate: they do not reproduce grid values.
Even with zero velocity, each step blurs the field with weights
(1/6, 4/6, 1/6) per dimension. This adds large numerical diffusion.
Fixing it means changing the interpolation in the advection scheme, so
the scheme replacement in Section 3 resolves it.

**Gravity has no effect.** `BodyForce` adds `-g * dt` to every J-face
velocity, but in a closed box of constant density a uniform force is the
gradient of a hydrostatic pressure, and the projection removes it to
solver tolerance (`tstVelocityCorrector`). The `-g` option of the
`advection` example therefore changes nothing. Gravity moves the fluid
only once the force depends on the advected quantity, for example through
a Boussinesq buoyancy term on the J-faces.

### Fixed bugs

The following bugs from the original review are fixed, each with a test
where one could be written:

- The pressure correction applied the J-face wall condition to `u`
  instead of `v` and converted face indices with the cell L2G map
  (`tstVelocityCorrector`).
- The RK3 back-trace computed its third stage from `v0`, so it was only
  second order; the 3D path also wrote the wrong stage (`tstTimeIntegrator`).
- Ghosts outside the physical domain were never filled. `BoundaryCondition::
  fillGhosts` now reflects each field across the walls after every gather
  (`tstBoundaryConditions`).
- The time step was fixed once from the inflow velocity, so back-traces
  could leave the halo. `TimeIntegrator::cflTimeStep` now limits each step
  to one cell of travel at the fastest face velocity, and the command-line
  time step is a maximum (`tstTimeIntegrator`).
- With all-solid walls the pressure system is singular. The corrector now
  removes the mean of the right-hand side (`tstVelocityCorrector`).
- The mesh divisibility check used an absolute tolerance (`tstMesh`).
- The inflow source tested only entity centers against its box, so a box
  smaller than a cell could miss every cell, and it forced a minimum
  quantity and velocity instead of adding a rate. It now weights each
  cell and face by the fraction of its control volume the box covers,
  adds `rate * dt` scaled by that fraction to the quantity, and blends
  face velocities toward the inflow velocity by the same fraction
  (`tstInflowSource`).
- The initializer received local indices; it now receives global ones
  (`tstProblemManager`).
- The halo test encoded local indices, the mesh test assumed even
  decompositions, and the boundary condition tests were empty; every test
  now runs at 1, 2, and 4 ranks.
- Untested fixes: unbalanced profiling regions, the HIP backend typo, the
  missing `bc` argument in the 3D advection call, the overwritten initial
  output and leaked `strdup` names, stale ghosts in the velocity output,
  command-line parsing and help text in `examples/advection.cpp`, and the
  README and example descriptions.

## 2. Unit Test Plan

### Infrastructure

`tests/CMakeLists.txt` hands each suite to Cabana's test harness in
`cmake/test_harness/`, which builds one executable per enabled Kokkos
backend and registers it at 1, 2, and 4 ranks through
`MPIEXEC_EXECUTABLE`; the GitHub workflow runs them all with `ctest`. A
new suite is a `tst<Name>.hpp` typed on `MeshDeviceTypes` from
`tstDriver.hpp`, added to `MPI_TESTS`; shared fixtures live in
`*Fixture.hpp`.

- Add a small helper header for tests: copy a `Cabana::Grid` array to the
  host, compute a global max, sum, or L2 norm over owned entities with
  `MPI_Allreduce`, and fill an array from a function of global position.
- Encode global, not local, indices in any value that crosses ranks.
- Keep test meshes small (32x32 to 128x128) so the suite finishes in
  seconds on CPU backends.

### Tests by component

The existing tests cover mesh decomposition and divisibility, halo
exchange, initializer indices, wall velocities, the pressure matrix, wall
ghost filling, RK3 order, the CFL time step, the pressure correction's
wall velocities, divergence, singular solve, and removal of uniform
gravity, and the inflow source's overlap weighting. The tests below
remain.

**ProblemManager**
- `advance` swaps current and next for every field, not only the cell
  field.
- The init functor receives the cell-center and face-center coordinates
  that `local_mesh.coordinates` reports.

**Interpolation**
- Linear interpolation reproduces any linear field exactly.
- The chosen higher-order interpolant reproduces grid values at grid
  points; a failing version of this test documents the B-spline bug in
  Section 1.
- Face-centered interpolation lands on the right staggered locations.

**TimeIntegrator**
- `rk3` with constant velocity returns `x0 - dt * u` exactly.
- `advect` with zero velocity leaves the field unchanged to round-off
  (fails until the B-spline bug is fixed).
- `advect` near a wall, at the largest time step `cflTimeStep` allows,
  stays in bounds and produces finite values. Run this under a Kokkos
  debug build with bounds checking.

**VelocityCorrector**
- Projection is idempotent: a second call changes nothing beyond the
  solver tolerance.
- Given `u = u_df + grad(phi)` for a known divergence-free `u_df`,
  projection recovers `u_df` to second order in `h`.

**Solver (integration tests)**
- Closed box, no inflow, no body force: total quantity is conserved to the
  scheme's tolerance, and the quantity stays within its initial min and
  max.
- One-rank and four-rank runs of the same short problem agree to round-off.
- Output: a run writes files for step 0 and every `write_freq` steps
  after it, and cell-centered velocities agree across rank counts.
- Command-line parsing: move `parseInput` into a testable function and
  check each option, including the validation fixed earlier.

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
- **Time step:** choose `dt = C * h / max|u|` with `C <= 0.8`. Add the
  factor `C` to `TimeIntegrator::cflTimeStep`, which already computes
  `max|u|` each step with a global reduction and uses `C = 1`.

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

1. **Land the remaining unit tests** from Section 2, so the new scheme
   starts from a tested base.
2. **Add flux-form advection** in a new `src/Advection.hpp`: a limiter
   function, a face-flux kernel for each direction, and a flux-divergence
   kernel. Document each function with the equation it implements and a
   textbook reference.
3. **Replace `TimeIntegrator::step`** with SSP-RK3. Store one extra stage
   array per field in `ProblemManager`, and gather and call
   `BoundaryCondition::fillGhosts` before each stage. `fillGhosts` gives
   free-slip walls; add a no-slip option that negates tangential velocity
   if viscosity is added.
4. **Reduce the halo to two cells** once nothing needs three, and update
   the comment in `ProblemManager.hpp` that explains the halo depth.
5. **Keep the old scheme for comparison** only if benchmark users need it.
   Otherwise remove it, along with the spline interpolation of advected
   quantities. Keep linear interpolation for output.
6. **Write `docs/numerics.md`:** the governing equations, the MAC grid
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
