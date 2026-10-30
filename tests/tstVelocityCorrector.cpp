#include <Cabana_Core.hpp>
#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>
#include <VelocityCorrector.hpp>

#include <mpi.h>

#include "tstDriver.hpp"
#include "tstProblemManager.hpp"

using Cell = Cabana::Grid::Cell;
using FaceI = Cabana::Grid::Face<Cabana::Grid::Dim::I>;
using FaceJ = Cabana::Grid::Face<Cabana::Grid::Dim::J>;
using Velocity = CabanaFluids::Field::Velocity;
using Current = CabanaFluids::Version::Current;

template <class T>
class VelocityCorrectorTest : public ::testing::Test
{
  protected:
    using ExecutionSpace = typename T::ExecutionSpace;
    using MemorySpace = typename T::MemorySpace;
    using mesh_type = CabanaFluids::Mesh<2, ExecutionSpace, MemorySpace>;
    using pm_type =
        CabanaFluids::ProblemManager<2, ExecutionSpace, MemorySpace>;

    const int boxCells_ = 64;

    void SetUp() override
    {
        Cabana::Grid::DimBlockPartitioner<2> partitioner;
        mesh_ = std::make_shared<mesh_type>(
            Kokkos::Array<double, 4>{ 0, 0, 1.0, 1.0 },
            std::array<int, 2>{ boxCells_, boxCells_ }, partitioner, 3,
            MPI_COMM_WORLD );
        pm_ = std::make_shared<pm_type>( mesh_, NullInitFunctor<2>() );
        bc_.boundary_type = { CabanaFluids::BoundaryType::SOLID,
                              CabanaFluids::BoundaryType::SOLID,
                              CabanaFluids::BoundaryType::SOLID,
                              CabanaFluids::BoundaryType::SOLID };
        bc_.min = mesh_->minDomainGlobalCellIndex();
        bc_.max = mesh_->maxDomainGlobalCellIndex();
    }

    // Set each owned I-face velocity to f(global i, global j)
    template <class Func_t>
    void setU( const Func_t& f )
    {
        auto local_grid = mesh_->localGrid();
        auto space = local_grid->indexSpace( Cabana::Grid::Own(), FaceI(),
                                             Cabana::Grid::Local() );
        auto l2g =
            Cabana::Grid::IndexConversion::createL2G( *local_grid, FaceI() );
        auto u = pm_->get( FaceI(), Velocity(), Current() );
        Kokkos::parallel_for(
            "set u", createExecutionPolicy( space, ExecutionSpace() ),
            KOKKOS_LAMBDA( const int i, const int j ) {
                int gi, gj;
                l2g( i, j, gi, gj );
                u( i, j, 0 ) = f( gi, gj );
            } );
    }

    void correct()
    {
        auto vc = CabanaFluids::createVelocityCorrector<2, ExecutionSpace,
                                                        MemorySpace>(
            pm_, bc_, 1.0, 0.01, "Reference", "Jacobi" );
        vc->correctVelocity();
    }

    std::shared_ptr<mesh_type> mesh_;
    std::shared_ptr<pm_type> pm_;
    CabanaFluids::BoundaryCondition<2> bc_;
};

TYPED_TEST_SUITE( VelocityCorrectorTest, MeshDeviceTypes );

// A single nonzero face in the middle of the domain: a divergence dipole
struct Dipole
{
    int n;
    KOKKOS_INLINE_FUNCTION double operator()( int gi, int gj ) const
    {
        return ( gi == n / 2 && gj == n / 2 ) ? 1.0 : 0.0;
    }
};

// Uniform unit divergence plus the dipole, so the right-hand side of the
// pressure solve has a nonzero mean
struct Expansion
{
    int n;
    KOKKOS_INLINE_FUNCTION double operator()( int gi, int gj ) const
    {
        return double( gi ) / n + Dipole{ n }( gi, gj );
    }
};

TYPED_TEST( VelocityCorrectorTest, SolidWallsAndDivergence )
{
    this->setU( Dipole{ this->boxCells_ } );
    this->correct();

    // Check that every solid-wall normal velocity is zero and that the
    // corrected field is divergence free
    auto pm = this->pm_;
    auto local_grid = this->mesh_->localGrid();
    pm->gather( Current() );
    auto u = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), pm->get( FaceI(), Velocity(), Current() ) );
    auto v = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), pm->get( FaceJ(), Velocity(), Current() ) );
    int n = this->boxCells_;

    auto checkWalls = [&]( auto entity, auto& view, int dim )
    {
        auto space = local_grid->indexSpace( Cabana::Grid::Own(), entity,
                                             Cabana::Grid::Local() );
        auto l2g =
            Cabana::Grid::IndexConversion::createL2G( *local_grid, entity );
        for ( int i = space.min( 0 ); i < space.max( 0 ); i++ )
            for ( int j = space.min( 1 ); j < space.max( 1 ); j++ )
            {
                int g[2];
                l2g( i, j, g[0], g[1] );
                if ( g[dim] == 0 || g[dim] == n )
                    ASSERT_EQ( view( i, j, 0 ), 0.0 )
                        << "wall face " << g[0] << ", " << g[1];
            }
    };
    checkWalls( FaceI(), u, 0 );
    checkWalls( FaceJ(), v, 1 );

    // The initial divergence peaks at n; the solver tolerance is 1e-6
    auto cells = local_grid->indexSpace( Cabana::Grid::Own(), Cell(),
                                         Cabana::Grid::Local() );
    for ( int i = cells.min( 0 ); i < cells.max( 0 ); i++ )
        for ( int j = cells.min( 1 ); j < cells.max( 1 ); j++ )
            ASSERT_LT( std::abs( u( i + 1, j, 0 ) - u( i, j, 0 ) +
                                 v( i, j + 1, 0 ) - v( i, j, 0 ) ) *
                           n,
                       1.0e-4 * n );
}

TYPED_TEST( VelocityCorrectorTest, SingularSystemConverges )
{
    // With every wall solid the pressure is defined only up to a constant,
    // and a right-hand side with nonzero mean has no solution unless the
    // corrector removes the mean
    this->setU( Expansion{ this->boxCells_ } );
    EXPECT_NO_THROW( this->correct() );
}
