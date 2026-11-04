#include <BodyForce.hpp>
#include <Cabana_Core.hpp>
#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>
#include <VelocityCorrector.hpp>

#include <mpi.h>

#include <string>

#include "ProblemManagerFixture.hpp"
#include "tstDriver.hpp"

using Cell = Cabana::Grid::Cell;
using FaceI = Cabana::Grid::Face<Cabana::Grid::Dim::I>;
using FaceJ = Cabana::Grid::Face<Cabana::Grid::Dim::J>;
using Velocity = CabanaFluids::Field::Velocity;
using Current = CabanaFluids::Version::Current;

// Each case pairs the device under test with a pressure solver and its
// preconditioner, named as createVelocityCorrector takes them
using TestDevice = DeviceType<TEST_EXECSPACE, TEST_MEMSPACE>;
struct ReferenceJacobi : TestDevice
{
    static constexpr const char* solver = "Reference";
    static constexpr const char* precon = "Jacobi";
};
#ifdef Cabana_ENABLE_HYPRE
struct PCGDiagonal : TestDevice
{
    static constexpr const char* solver = "PCG";
    static constexpr const char* precon = "Diagonal";
};
struct PCGPFMG : TestDevice
{
    static constexpr const char* solver = "PCG";
    static constexpr const char* precon = "PFMG";
};
struct GMRESPFMG : TestDevice
{
    static constexpr const char* solver = "GMRES";
    static constexpr const char* precon = "PFMG";
};
struct BiCGSTABDiagonal : TestDevice
{
    static constexpr const char* solver = "BiCGSTAB";
    static constexpr const char* precon = "Diagonal";
};
struct PFMGNone : TestDevice
{
    static constexpr const char* solver = "PFMG";
    static constexpr const char* precon = "none";
};
struct SMGNone : TestDevice
{
    static constexpr const char* solver = "SMG";
    static constexpr const char* precon = "none";
};
#endif

#ifdef Cabana_ENABLE_HYPRE
using SolverCases =
    ::testing::Types<ReferenceJacobi, PCGDiagonal, PCGPFMG, GMRESPFMG,
                     BiCGSTABDiagonal, PFMGNone, SMGNone>;
#else
using SolverCases = ::testing::Types<ReferenceJacobi>;
#endif

struct SolverCaseName
{
    template <class T>
    static std::string GetName( int )
    {
        return std::string( T::solver ) + T::precon;
    }
};

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
        setU( pm_, f );
    }

    template <class Func_t>
    void setU( const std::shared_ptr<pm_type>& pm, const Func_t& f )
    {
        auto local_grid = mesh_->localGrid();
        auto space = local_grid->indexSpace( Cabana::Grid::Own(), FaceI(),
                                             Cabana::Grid::Local() );
        auto l2g =
            Cabana::Grid::IndexConversion::createL2G( *local_grid, FaceI() );
        auto u = pm->get( FaceI(), Velocity(), Current() );
        Kokkos::parallel_for(
            "set u", createExecutionPolicy( space, ExecutionSpace() ),
            KOKKOS_LAMBDA( const int i, const int j ) {
                int gi, gj;
                l2g( i, j, gi, gj );
                u( i, j, 0 ) = f( gi, gj );
            } );
    }

    void correct() { correct( pm_, T::solver, T::precon ); }

    void correct( const std::shared_ptr<pm_type>& pm, const std::string& solver,
                  const std::string& precon )
    {
        auto vc = CabanaFluids::createVelocityCorrector<2, ExecutionSpace,
                                                        MemorySpace>(
            pm, bc_, 1.0, 0.01, solver, precon );
        vc->correctVelocity();
    }

    std::shared_ptr<mesh_type> mesh_;
    std::shared_ptr<pm_type> pm_;
    CabanaFluids::BoundaryCondition<2> bc_;
};

TYPED_TEST_SUITE( VelocityCorrectorTest, SolverCases, SolverCaseName );

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

TYPED_TEST( VelocityCorrectorTest, UniformGravityProjectsAway )
{
    // Uniform gravity on a constant-density fluid in a closed box is the
    // gradient of a hydrostatic pressure, so the projection removes it and
    // gravity alone never moves the fluid
    using ExecutionSpace = typename TypeParam::ExecutionSpace;
    auto local_grid = this->mesh_->localGrid();
    auto space = local_grid->indexSpace( Cabana::Grid::Own(), FaceJ(),
                                         Cabana::Grid::Local() );
    auto l2g = Cabana::Grid::IndexConversion::createL2G( *local_grid, FaceJ() );
    auto v = this->pm_->get( FaceJ(), Velocity(), Current() );
    CabanaFluids::BodyForce<2> gravity( 0.0, -9.8 );
    auto bc = this->bc_;
    Kokkos::parallel_for(
        "apply gravity", createExecutionPolicy( space, ExecutionSpace() ),
        KOKKOS_LAMBDA( const int i, const int j ) {
            int gi, gj;
            l2g( i, j, gi, gj );
            gravity( FaceJ(), v, i, j, 0.0, 0.0, 0.01, 0.0 );
            bc( FaceJ(), v, gi, gj, i, j );
        } );
    this->correct();

    auto v_host = Kokkos::create_mirror_view_and_copy( Kokkos::HostSpace(), v );
    for ( int i = space.min( 0 ); i < space.max( 0 ); i++ )
        for ( int j = space.min( 1 ); j < space.max( 1 ); j++ )
            ASSERT_LT( std::abs( v_host( i, j, 0 ) ), 1.0e-4 * 9.8 * 0.01 );
}

TYPED_TEST( VelocityCorrectorTest, MatchesReferenceSolve )
{
    if ( std::string( TypeParam::solver ) == "Reference" )
        GTEST_SKIP() << "compares each hypre solver with the reference";

    // Project the same dipole with the reference solver and with this
    // case's solver. The pressures may differ by a constant, since the
    // all-solid system is singular, but the corrected velocities may not.
    using pm_type = typename TestFixture::pm_type;
    auto ref_pm =
        std::make_shared<pm_type>( this->mesh_, NullInitFunctor<2>() );
    this->setU( Dipole{ this->boxCells_ } );
    this->setU( ref_pm, Dipole{ this->boxCells_ } );
    this->correct();
    this->correct( ref_pm, "Reference", "Jacobi" );

    auto local_grid = this->mesh_->localGrid();
    auto compare = [&]( auto entity )
    {
        auto x = Kokkos::create_mirror_view_and_copy(
            Kokkos::HostSpace(),
            this->pm_->get( entity, Velocity(), Current() ) );
        auto ref = Kokkos::create_mirror_view_and_copy(
            Kokkos::HostSpace(), ref_pm->get( entity, Velocity(), Current() ) );
        auto space = local_grid->indexSpace( Cabana::Grid::Own(), entity,
                                             Cabana::Grid::Local() );
        // The initial velocity peaks at 1 and both solvers stop at a
        // relative residual of 1e-6, so allow ten times that
        for ( int i = space.min( 0 ); i < space.max( 0 ); i++ )
            for ( int j = space.min( 1 ); j < space.max( 1 ); j++ )
                ASSERT_NEAR( x( i, j, 0 ), ref( i, j, 0 ), 1.0e-5 )
                    << "local face " << i << ", " << j;
    };
    compare( FaceI() );
    compare( FaceJ() );
}
