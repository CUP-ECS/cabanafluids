#include <Cabana_Core.hpp>
#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>
#include <TimeIntegrator.hpp>

#include <mpi.h>

#include "tstDriver.hpp"
#include "tstMesh.hpp"
#include "tstProblemManager.hpp"

TYPED_TEST_SUITE( ProblemManagerTest, MeshDeviceTypes );

using FaceI = Cabana::Grid::Face<Cabana::Grid::Dim::I>;
using FaceJ = Cabana::Grid::Face<Cabana::Grid::Dim::J>;
using Velocity = CabanaFluids::Field::Velocity;
using Current = CabanaFluids::Version::Current;

// Fill a face velocity component over the ghosted space with a function of
// position
template <class ExecutionSpace, class Entity_t, class PM_t, class Func_t>
void fillVelocity( const PM_t& pm, Entity_t entity, const Func_t& func )
{
    auto local_grid = pm->mesh()->localGrid();
    auto local_mesh = *( pm->mesh()->localMesh() );
    auto u = pm->get( entity, Velocity(), Current() );
    auto space = local_grid->indexSpace( Cabana::Grid::Ghost(), entity,
                                         Cabana::Grid::Local() );
    Kokkos::parallel_for(
        "fill velocity", createExecutionPolicy( space, ExecutionSpace() ),
        KOKKOS_LAMBDA( const int i, const int j ) {
            int idx[2] = { i, j };
            double x[2];
            local_mesh.coordinates( entity, idx, x );
            u( i, j, 0 ) = func( x );
        } );
}

struct RotationU
{
    double c[2];
    KOKKOS_INLINE_FUNCTION double operator()( const double x[2] ) const
    {
        return -( x[1] - c[1] );
    }
};

struct RotationV
{
    double c[2];
    KOKKOS_INLINE_FUNCTION double operator()( const double x[2] ) const
    {
        return x[0] - c[0];
    }
};

TYPED_TEST( ProblemManagerTest, RK3ThirdOrder )
{
    using ExecutionSpace = typename TestFixture::ExecutionSpace;

    // Rotate counterclockwise at unit angular velocity about the center of
    // this rank's owned cells. Linear interpolation reproduces this linear
    // field exactly, so the trace error is the integrator's error alone.
    auto pm = this->testPM_;
    auto mesh = pm->mesh();
    auto& global_grid = mesh->localGrid()->globalGrid();
    double h = mesh->cellSize();
    double c[2];
    for ( int d = 0; d < 2; d++ )
        c[d] = ( global_grid.globalOffset( d ) +
                 0.5 * global_grid.ownedNumCell( d ) ) *
               h;
    fillVelocity<ExecutionSpace>( pm, FaceI(), RotationU{ { c[0], c[1] } } );
    fillVelocity<ExecutionSpace>( pm, FaceJ(), RotationV{ { c[0], c[1] } } );

    auto u = pm->get( FaceI(), Velocity(), Current() );
    auto v = pm->get( FaceJ(), Velocity(), Current() );
    auto local_mesh = *( mesh->localMesh() );
    double r = 2 * h;
    auto traceError = [&]( double dt )
    {
        double error;
        Kokkos::parallel_reduce(
            "rk3 trace", Kokkos::RangePolicy<ExecutionSpace>( 0, 1 ),
            KOKKOS_LAMBDA( const int, double& e ) {
                double x0[2] = { c[0] + r, c[1] };
                double trace[2];
                CabanaFluids::TimeIntegrator::rk3<2>( x0, local_mesh, u, v, dt,
                                                      trace );
                e = Kokkos::hypot( trace[0] - ( c[0] + r * Kokkos::cos( dt ) ),
                                   trace[1] -
                                       ( c[1] - r * Kokkos::sin( dt ) ) );
            },
            error );
        return error;
    };

    // A third-order method has fourth-order local error, so halving the
    // step cuts the error by about 16; a second-order one by about 8
    EXPECT_GT( traceError( 0.5 ) / traceError( 0.25 ), 12.0 );
}

TYPED_TEST( ProblemManagerTest, CFLTimeStep )
{
    using ExecutionSpace = typename TestFixture::ExecutionSpace;

    auto pm = this->testPM_;
    auto mesh = pm->mesh();
    double h = mesh->cellSize();
    auto u = pm->get( FaceI(), Velocity(), Current() );
    auto v = pm->get( FaceJ(), Velocity(), Current() );
    Kokkos::deep_copy( u, 0.0 );
    Kokkos::deep_copy( v, 0.0 );

    // No flow leaves the requested maximum step
    EXPECT_EQ(
        CabanaFluids::TimeIntegrator::cflTimeStep( ExecutionSpace(), *pm, 1.0 ),
        1.0 );

    // The fastest face on any rank limits the step on every rank
    if ( mesh->rank() == 0 )
    {
        auto space = mesh->localGrid()->indexSpace(
            Cabana::Grid::Own(), FaceJ(), Cabana::Grid::Local() );
        Kokkos::deep_copy(
            Kokkos::subview( v, space.min( 0 ), space.min( 1 ), 0 ), -4.0 );
    }
    EXPECT_DOUBLE_EQ(
        CabanaFluids::TimeIntegrator::cflTimeStep( ExecutionSpace(), *pm, 1.0 ),
        h / 4.0 );
    EXPECT_EQ( CabanaFluids::TimeIntegrator::cflTimeStep( ExecutionSpace(), *pm,
                                                          h / 8.0 ),
               h / 8.0 );
}
