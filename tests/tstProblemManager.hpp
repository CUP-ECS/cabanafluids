#include <Cabana_Core.hpp>
#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>
#include <ProblemManager.hpp>

#include <mpi.h>

#include "MeshFixture.hpp"
#include "ProblemManagerFixture.hpp"
#include "tstDriver.hpp"

TYPED_TEST_SUITE( ProblemManagerTest, MeshDeviceTypes );

using Cell = Cabana::Grid::Cell;
using FaceI = Cabana::Grid::Face<Cabana::Grid::Dim::I>;
using FaceJ = Cabana::Grid::Face<Cabana::Grid::Dim::J>;
using FaceK = Cabana::Grid::Face<Cabana::Grid::Dim::K>;
using Quantity = CabanaFluids::Field::Quantity;
using Velocity = CabanaFluids::Field::Velocity;
using Current = CabanaFluids::Version::Current;
using Next = CabanaFluids::Version::Next;

TYPED_TEST( ProblemManagerTest, StateArrayTest )
{
    using ExecutionSpace = typename TestFixture::ExecutionSpace;

    // Get basic mesh state
    auto pm = this->testPM_;
    auto mesh = pm->mesh();
    auto qcurr = pm->get( Cell(), Quantity(), Current() );
    auto qnext = pm->get( Cell(), Quantity(), Next() );
    auto rank = mesh->rank();

    /* Set values in the array based on our rank. Each cell gets a value of
     * rank*1000 + i * 100 + j * 10 + an int based on its type as follows
     * FaceI: +1
     * FaceJ: +2
     * Next(): +5
     */
    auto qspace = mesh->localGrid()->indexSpace( Cabana::Grid::Own(), Cell(),
                                                 Cabana::Grid::Local() );
    Kokkos::parallel_for(
        "InitializeCellFields",
        createExecutionPolicy( qspace, ExecutionSpace() ),
        KOKKOS_LAMBDA( const int i, const int j ) {
            qcurr( i, j, 0 ) = rank * 1000 + i * 100 + j * 10;
            qnext( i, j, 0 ) = rank * 1000 + i * 100 + j * 10 + 5;
        } );

    // Check that we can swap the views properly by checking the cell view
    // (we don't check the other views for now)
    pm->advance( Cell(), Quantity() );
    qcurr = pm->get( Cell(), Quantity(), Current() );
    auto qcopy =
        Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), qcurr );
    for ( int i = qspace.min( 0 ); i < qspace.max( 0 ); i++ )
        for ( int j = qspace.min( 1 ); j < qspace.max( 1 ); j++ )
            ASSERT_EQ( qcopy( i, j, 0 ), rank * 1000 + i * 100 + j * 10 + 5 );
}

TYPED_TEST( ProblemManagerTest, HaloTest )
{
    using ExecutionSpace = typename TestFixture::ExecutionSpace;

    auto pm = this->testPM_;
    auto mesh = pm->mesh();
    auto ucurr = pm->get( FaceI(), Velocity(), Current() );
    auto uspace = mesh->localGrid()->indexSpace( Cabana::Grid::Own(), FaceI(),
                                                 Cabana::Grid::Local() );
    auto l2g =
        Cabana::Grid::IndexConversion::createL2G( *mesh->localGrid(), FaceI() );

    // Encode global indices, which every rank agrees on, in the values
    Kokkos::parallel_for(
        "InitializeFaceIFields",
        createExecutionPolicy( uspace, ExecutionSpace() ),
        KOKKOS_LAMBDA( const int i, const int j ) {
            int gi, gj;
            l2g( i, j, gi, gj );
            ucurr( i, j, 0 ) = gi * 10000 + gj;
        } );

    // Check that we can halo the views appropriately, using the FaceI direction
    // as the check
    pm->gather( Current() );
    std::array<int, 2> directions[4] = {
        { -1, 0 }, { 1, 0 }, { 0, -1 }, { 0, 1 } };
    auto ucopy =
        Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), ucurr );
    for ( int d = 0; d < 4; d++ )
    {
        auto dir = directions[d];
        if ( mesh->localGrid()->neighborRank( dir ) < 0 )
            continue;
        auto u_shared_space = mesh->localGrid()->sharedIndexSpace(
            Cabana::Grid::Ghost(), FaceI(), dir );
        for ( int i = u_shared_space.min( 0 ); i < u_shared_space.max( 0 );
              i++ )
            for ( int j = u_shared_space.min( 1 ); j < u_shared_space.max( 1 );
                  j++ )
            {
                int gi, gj;
                l2g( i, j, gi, gj );
                ASSERT_EQ( ucopy( i, j, 0 ), gi * 10000 + gj );
            }
    }
}

// Records the index each entity receives from the initializer
template <std::size_t Dim>
class IndexInitFunctor
{
  public:
    template <class Entity_t, class Field_t>
    KOKKOS_INLINE_FUNCTION bool
    operator()( Entity_t, Field_t, const int index[Dim],
                [[maybe_unused]] const double x[Dim], double& value ) const
    {
        value = index[0] * 10000 + index[1];
        return true;
    };
};

template <class Entity_t, class Field_t, class PM_t>
void checkInitIndices( const PM_t& pm, Entity_t entity, Field_t field )
{
    auto local_grid = pm->mesh()->localGrid();
    auto space = local_grid->indexSpace( Cabana::Grid::Own(), entity,
                                         Cabana::Grid::Local() );
    auto l2g = Cabana::Grid::IndexConversion::createL2G( *local_grid, entity );
    auto values = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), pm->get( entity, field, Current() ) );
    for ( int i = space.min( 0 ); i < space.max( 0 ); i++ )
        for ( int j = space.min( 1 ); j < space.max( 1 ); j++ )
        {
            int gi, gj;
            l2g( i, j, gi, gj );
            ASSERT_EQ( values( i, j, 0 ), gi * 10000 + gj );
        }
}

TYPED_TEST( ProblemManagerTest, InitializerGetsGlobalIndices )
{
    using pm_type =
        CabanaFluids::ProblemManager<2, typename TestFixture::ExecutionSpace,
                                     typename TestFixture::MemorySpace>;
    auto pm =
        std::make_shared<pm_type>( this->testMesh_, IndexInitFunctor<2>() );

    checkInitIndices( pm, Cell(), Quantity() );
    checkInitIndices( pm, FaceI(), Velocity() );
    checkInitIndices( pm, FaceJ(), Velocity() );
}
