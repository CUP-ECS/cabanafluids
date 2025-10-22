# We use the cabana test harness here for test cases. It is in root/cmake/test_harness

set(MPI_TESTS
  Mesh
  ProblemManager
  BoundaryConditions
  )

Cabana_add_tests(MPI PACKAGE cabanafluids NAMES ${MPI_TESTS})
