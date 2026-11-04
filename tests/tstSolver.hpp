#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>
#include <Solver.hpp>

#include <mpi.h>

#include "ProblemManagerFixture.hpp"
#include "tstDriver.hpp"

// The example selects its Kokkos backend by name, so a build that enables
// Threads must accept "threads" in createSolver
#ifdef KOKKOS_ENABLE_THREADS
TEST( Solver, CreatesThreadsSolver )
{
    Cabana::Grid::DimBlockPartitioner<2> partitioner;
    CabanaFluids::BoundaryCondition<2> bc;
    bc.boundary_type = {
        CabanaFluids::BoundaryType::SOLID, CabanaFluids::BoundaryType::SOLID,
        CabanaFluids::BoundaryType::SOLID, CabanaFluids::BoundaryType::SOLID };
    EXPECT_NO_THROW( CabanaFluids::createSolver(
        "threads", MPI_COMM_WORLD, { 0, 0, 1.0, 1.0 }, { 16, 16 }, partitioner,
        1.0, NullInitFunctor<2>(), bc,
        CabanaFluids::InflowSource<2>( { 0.5, 0.5 }, { 0.1, 0.1 }, { 0.0, 0.0 },
                                       0.0 ),
        CabanaFluids::BodyForce<2>( 0.0, 0.0 ), 0.01, "Reference", "Jacobi" ) );
}
#endif
