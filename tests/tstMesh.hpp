#include "gtest/gtest.h"

#include <Cabana_Core.hpp>
#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <Mesh.hpp>

#include <mpi.h>

#include "MeshFixture.hpp"
#include "tstDriver.hpp"

TYPED_TEST_SUITE( MeshTest, MeshDeviceTypes );

TYPED_TEST( MeshTest, BasicParameters )
{
    int r;
    EXPECT_EQ( this->testMesh_->cellSize(), this->boxWidth_ / this->boxCells_ );

    auto mins = this->testMesh_->minDomainGlobalCellIndex();
    EXPECT_EQ( mins[0], 0 );
    EXPECT_EQ( mins[1], 0 );
    auto maxs = this->testMesh_->maxDomainGlobalCellIndex();
    EXPECT_EQ( maxs[0], this->boxCells_ - 1 );
    EXPECT_EQ( maxs[1], this->boxCells_ - 1 );

    MPI_Comm_rank( MPI_COMM_WORLD, &r );
    EXPECT_EQ( this->testMesh_->rank(), r );
}

TYPED_TEST( MeshTest, LocalGridSetup )
{
    /* Here we check that the local grid is decomposed like
     * we think it should be. That is, the number of ghosts cells
     * is right, the index spaces for owned, ghost, and boundary
     * cells are right, and so on. */
    auto local_grid = this->testMesh_->localGrid();
    auto global_grid = local_grid->globalGrid();
    for ( int i = 0; i < 2; i++ )
    {
        EXPECT_EQ( this->boxCells_,
                   global_grid.globalNumEntity( Cabana::Grid::Cell(), i ) );
    }

    /* Make sure the number of owned cells is our share of what was requested */
    auto own_local_cell_space = local_grid->indexSpace(
        Cabana::Grid::Own(), Cabana::Grid::Cell(), Cabana::Grid::Local() );
    for ( int i = 0; i < 2; i++ )
    {
        int share = this->boxCells_ / global_grid.dimNumBlock( i );
        EXPECT_GE( own_local_cell_space.extent( i ), share );
        EXPECT_LE( own_local_cell_space.extent( i ), share + 1 );
    }

    /*
     * Next we extract the ghosted faces, which encompass the owned faces and
     * the ghosts in each dimension. Note that for faces there is a template
     * dimension parameter - there is a * separate index space for each
     * spatial dimension.
     */
    auto ghost_local_face_space = local_grid->indexSpace(
        Cabana::Grid::Ghost(), Cabana::Grid::Face<Cabana::Grid::Dim::I>(),
        Cabana::Grid::Local() );
    EXPECT_EQ( ghost_local_face_space.extent( 0 ),
               own_local_cell_space.extent( 0 ) + 2 * this->haloWidth_ + 1 );
    EXPECT_EQ( ghost_local_face_space.extent( 1 ),
               own_local_cell_space.extent( 1 ) + 2 * this->haloWidth_ );
}

TYPED_TEST( MeshTest, DivisibilityCheck )
{
    using mesh_type = CabanaFluids::Mesh<2, typename TypeParam::ExecutionSpace,
                                         typename TypeParam::MemorySpace>;

    // 88 cells of 100 / 88 miss the extent by 1.4e-14: more than ten
    // epsilons, but well within the relative round-off of the extent
    Kokkos::Array<double, 4> large_box = { 0, 0, 100.0, 100.0 };
    EXPECT_NO_THROW( mesh_type( large_box, { 88, 88 }, this->partitioner_,
                                this->haloWidth_, MPI_COMM_WORLD ) );

    Kokkos::Array<double, 4> uneven_box = { 0, 0, 1.0, 1.5 };
    EXPECT_THROW( mesh_type( uneven_box, { 32, 32 }, this->partitioner_,
                             this->haloWidth_, MPI_COMM_WORLD ),
                  std::logic_error );
}
