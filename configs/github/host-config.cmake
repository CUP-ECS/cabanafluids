# OpenMPI needs an extra argument to run as root and needs to be able
# to oversubscribe CPUs to run multi-rank tests
set(MPIEXEC_PREFLAGS "--allow-run-as-root;--oversubscribe" CACHE STRING "")
