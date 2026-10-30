#include "gtest/gtest.h"

// Include Statements

#include <Cabana_Core.hpp>
#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <BoundaryConditions.hpp>

#include "tstBoundaryConditions.hpp"
#include "tstDriver.hpp"
#include "tstMesh.hpp"
#include "tstProblemManager.hpp"

TYPED_TEST_SUITE( BoundaryConditionsTest, MeshDeviceTypes );

using Cell = Cabana::Grid::Cell;
using FaceI = Cabana::Grid::Face<Cabana::Grid::Dim::I>;
using FaceJ = Cabana::Grid::Face<Cabana::Grid::Dim::J>;
using Quantity = CabanaFluids::Field::Quantity;
using Velocity = CabanaFluids::Field::Velocity;
using Current = CabanaFluids::Version::Current;
using CabanaFluids::BoundaryType;

template <class Mesh_t>
CabanaFluids::BoundaryCondition<2> makeBC( const Mesh_t& mesh, int type )
{
    CabanaFluids::BoundaryCondition<2> bc;
    bc.boundary_type = { type, type, type, type };
    bc.min = mesh->minDomainGlobalCellIndex();
    bc.max = mesh->maxDomainGlobalCellIndex();
    return bc;
}

// Set every owned face velocity to 1, apply the boundary condition, and
// return the result on the host
template <class ExecutionSpace, class Entity_t, class PM_t>
auto applyVelocityBC( const PM_t& pm,
                      const CabanaFluids::BoundaryCondition<2>& bc,
                      Entity_t entity )
{
    auto local_grid = pm->mesh()->localGrid();
    auto space = local_grid->indexSpace( Cabana::Grid::Own(), entity,
                                         Cabana::Grid::Local() );
    auto l2g = Cabana::Grid::IndexConversion::createL2G( *local_grid, entity );
    auto u = pm->get( entity, Velocity(), Current() );
    Kokkos::parallel_for(
        "apply velocity bc", createExecutionPolicy( space, ExecutionSpace() ),
        KOKKOS_LAMBDA( const int i, const int j ) {
            int gi, gj;
            l2g( i, j, gi, gj );
            u( i, j, 0 ) = 1.0;
            bc( entity, u, gi, gj, i, j );
        } );
    return Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), u );
}

// Check that only wall-normal velocities on solid walls are zeroed
template <class ExecutionSpace, class Entity_t, class PM_t>
void checkVelocityBC( const PM_t& pm, int type, Entity_t entity, int dim )
{
    auto bc = makeBC( pm->mesh(), type );
    auto u = applyVelocityBC<ExecutionSpace>( pm, bc, entity );
    auto local_grid = pm->mesh()->localGrid();
    auto space = local_grid->indexSpace( Cabana::Grid::Own(), entity,
                                         Cabana::Grid::Local() );
    auto l2g = Cabana::Grid::IndexConversion::createL2G( *local_grid, entity );
    for ( int i = space.min( 0 ); i < space.max( 0 ); i++ )
        for ( int j = space.min( 1 ); j < space.max( 1 ); j++ )
        {
            int g[2];
            l2g( i, j, g[0], g[1] );
            bool wall = g[dim] == bc.min[dim] || g[dim] == bc.max[dim] + 1;
            bool zero = wall && type == BoundaryType::SOLID;
            ASSERT_EQ( u( i, j, 0 ), zero ? 0.0 : 1.0 );
        }
}

TYPED_TEST( BoundaryConditionsTest, MeshSolidEdge2D )
{
    using ExecutionSpace = typename TestFixture::ExecutionSpace;
    checkVelocityBC<ExecutionSpace>( this->testPM_, BoundaryType::SOLID,
                                     FaceI(), 0 );
    checkVelocityBC<ExecutionSpace>( this->testPM_, BoundaryType::SOLID,
                                     FaceJ(), 1 );
}

TYPED_TEST( BoundaryConditionsTest, MeshFreeEdge2D )
{
    using ExecutionSpace = typename TestFixture::ExecutionSpace;
    checkVelocityBC<ExecutionSpace>( this->testPM_, BoundaryType::FREE, FaceI(),
                                     0 );
    checkVelocityBC<ExecutionSpace>( this->testPM_, BoundaryType::FREE, FaceJ(),
                                     1 );
}

// Build the 5-point pressure matrix with unit scale and check each entry
// against the expected weights for interior, edge, and corner cells
template <class ExecutionSpace, class MemorySpace, class Mesh_t>
void checkMatrix( const Mesh_t& mesh, int type )
{
    auto bc = makeBC( mesh, type );
    auto local_grid = mesh->localGrid();
    auto layout = Cabana::Grid::createArrayLayout( local_grid, 5, Cell() );
    auto entries =
        Cabana::Grid::createArray<double, MemorySpace>( "matrix", layout );
    auto m = entries->view();
    auto space = local_grid->indexSpace( Cabana::Grid::Own(), Cell(),
                                         Cabana::Grid::Local() );
    auto l2g = Cabana::Grid::IndexConversion::createL2G( *local_grid, Cell() );
    Kokkos::parallel_for(
        "fill matrix", createExecutionPolicy( space, ExecutionSpace() ),
        KOKKOS_LAMBDA( const int i, const int j ) {
            int gi, gj;
            l2g( i, j, gi, gj );
            m( i, j, 0 ) = 4.0;
            for ( int n = 1; n < 5; n++ )
                m( i, j, n ) = -1.0;
            bc.build_matrix( gi, gj, i, j, m, 1.0 );
        } );
    auto host = Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), m );

    for ( int i = space.min( 0 ); i < space.max( 0 ); i++ )
        for ( int j = space.min( 1 ); j < space.max( 1 ); j++ )
        {
            int gi, gj;
            l2g( i, j, gi, gj );
            // Stencil order: center, -i, +i, -j, +j
            bool outside[5] = { false, gi == bc.min[0], gi == bc.max[0],
                                gj == bc.min[1], gj == bc.max[1] };
            double diag = 4.0;
            for ( int n = 1; n < 5; n++ )
            {
                ASSERT_EQ( host( i, j, n ), outside[n] ? 0.0 : -1.0 );
                if ( outside[n] && type == BoundaryType::SOLID )
                    diag -= 1.0;
            }
            ASSERT_EQ( host( i, j, 0 ), diag );
        }
}

TYPED_TEST( BoundaryConditionsTest, MatrixSolidEdge2D )
{
    checkMatrix<typename TestFixture::ExecutionSpace,
                typename TestFixture::MemorySpace>( this->testMesh_,
                                                    BoundaryType::SOLID );
}

TYPED_TEST( BoundaryConditionsTest, MatrixFreeEdge2D )
{
    checkMatrix<typename TestFixture::ExecutionSpace,
                typename TestFixture::MemorySpace>( this->testMesh_,
                                                    BoundaryType::FREE );
}

// Write a nonzero function of global index into owned entities, gather and
// fill ghosts, and check that every ghost outside the domain mirrors the
// owned entity across the wall. Wall-normal velocity at a solid wall flips
// sign so that it interpolates to zero on the wall.
template <class ExecutionSpace, class Entity_t, class Field_t, class PM_t>
void checkGhostFill( const PM_t& pm, int type, Entity_t entity, Field_t field,
                     const bool face_aligned[2] )
{
    auto bc = makeBC( pm->mesh(), type );
    auto local_grid = pm->mesh()->localGrid();
    auto l2g = Cabana::Grid::IndexConversion::createL2G( *local_grid, entity );
    auto f = pm->get( entity, field, Current() );
    auto owned = local_grid->indexSpace( Cabana::Grid::Own(), entity,
                                         Cabana::Grid::Local() );
    Kokkos::parallel_for(
        "fill owned", createExecutionPolicy( owned, ExecutionSpace() ),
        KOKKOS_LAMBDA( const int i, const int j ) {
            int gi, gj;
            l2g( i, j, gi, gj );
            f( i, j, 0 ) = ( gi + 1000 ) * 10000 + gj + 1000;
        } );
    pm->gather( Current() );
    bc.fillGhosts( ExecutionSpace(), *local_grid, entity, f );

    auto host = Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), f );
    auto ghosted = local_grid->indexSpace( Cabana::Grid::Ghost(), entity,
                                           Cabana::Grid::Local() );
    for ( int i = ghosted.min( 0 ); i < ghosted.max( 0 ); i++ )
        for ( int j = ghosted.min( 1 ); j < ghosted.max( 1 ); j++ )
        {
            int g[2];
            l2g( i, j, g[0], g[1] );
            int src[2] = { g[0], g[1] };
            double sign = 1.0;
            bool outside = false;
            for ( int d = 0; d < 2; d++ )
            {
                int hi = bc.max[d] + ( face_aligned[d] ? 1 : 0 );
                int lo_mirror = face_aligned[d] ? 0 : -1;
                int hi_mirror = 2 * hi + ( face_aligned[d] ? 0 : 1 );
                if ( g[d] < bc.min[d] )
                    src[d] = lo_mirror - g[d];
                else if ( g[d] > hi )
                    src[d] = hi_mirror - g[d];
                else
                    continue;
                outside = true;
                if ( face_aligned[d] && type == BoundaryType::SOLID )
                    sign = -sign;
            }
            if ( outside )
                ASSERT_EQ( host( i, j, 0 ), sign * ( ( src[0] + 1000 ) * 10000 +
                                                     src[1] + 1000 ) )
                    << "global index " << g[0] << ", " << g[1];
        }
}

TYPED_TEST( BoundaryConditionsTest, GhostFill2D )
{
    using ExecutionSpace = typename TestFixture::ExecutionSpace;
    const bool cell[2] = { false, false };
    const bool face_i[2] = { true, false };
    const bool face_j[2] = { false, true };
    for ( int type : { BoundaryType::SOLID, BoundaryType::FREE } )
    {
        checkGhostFill<ExecutionSpace>( this->testPM_, type, Cell(), Quantity(),
                                        cell );
        checkGhostFill<ExecutionSpace>( this->testPM_, type, FaceI(),
                                        Velocity(), face_i );
        checkGhostFill<ExecutionSpace>( this->testPM_, type, FaceJ(),
                                        Velocity(), face_j );
    }
}
