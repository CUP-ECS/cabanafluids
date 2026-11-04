# OpenMPI needs an extra argument to run as root and needs to be able
# to oversubscribe CPUs to run multi-rank tests
set(MPIEXEC_PREFLAGS "--allow-run-as-root;--oversubscribe" CACHE STRING "")

# The container builds Cabana with hypre; fail rather than silently compile
# the hypre tests out if it ever stops doing so
set(CabanaFluids_REQUIRE_HYPRE ON CACHE BOOL "")
