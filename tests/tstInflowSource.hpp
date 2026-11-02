#include <Cabana_Grid.hpp>
#include <InflowSource.hpp>
#include <Kokkos_Core.hpp>

#include <mpi.h>

#include "tstDriver.hpp"

using Cell = Cabana::Grid::Cell;
using FaceI = Cabana::Grid::Face<Cabana::Grid::Dim::I>;
using FaceJ = Cabana::Grid::Face<Cabana::Grid::Dim::J>;
using HostArray = Kokkos::View<double***, Kokkos::HostSpace>;

// Apply the source on the unit square with n x n cells of size h, setting
// every entity to `initial` first. Entity centers follow the MAC layout:
// cells at ((i+1/2)h, (j+1/2)h), I-faces at (ih, (j+1/2)h), and J-faces at
// ((i+1/2)h, jh).
constexpr int n = 10;
constexpr double h = 1.0 / n;
constexpr double dt = 0.01;

template <class Entity>
HostArray apply( const CabanaFluids::InflowSource<2>& source, Entity,
                 double xoffset, double yoffset, double initial )
{
    HostArray a( "a", n + 1, n + 1, 1 );
    Kokkos::deep_copy( a, initial );
    for ( int i = 0; i <= n; i++ )
        for ( int j = 0; j <= n; j++ )
            source( Entity(), a, i, j, ( i + xoffset ) * h, ( j + yoffset ) * h,
                    dt, h * h );
    return a;
}

double injected( const HostArray& q )
{
    double total = 0.0;
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < n; j++ )
            total += q( i, j, 0 ) * h * h;
    return total;
}

// A box smaller than a cell that contains no cell center still deposits
// all of its inflow, in the cell that holds it
TEST( InflowSource, TinyBoxHitsACell )
{
    double rate = 3.0;
    CabanaFluids::InflowSource<2> source( { 0.42, 0.42 }, { 0.01, 0.01 },
                                          { 0.0, 0.0 }, rate );
    auto q = apply( source, Cell(), 0.5, 0.5, 0.0 );
    EXPECT_GT( q( 4, 4, 0 ), 0.0 );
    EXPECT_NEAR( injected( q ), rate * dt * 0.01 * 0.01, 1.0e-15 );
}

// A box that cuts through cells deposits rate * dt * area in total,
// split by how much of each cell the box covers
TEST( InflowSource, InflowIsConserved )
{
    double rate = 3.0;
    CabanaFluids::InflowSource<2> source( { 0.23, 0.37 }, { 0.31, 0.25 },
                                          { 0.0, 0.0 }, rate );
    auto q = apply( source, Cell(), 0.5, 0.5, 0.0 );
    EXPECT_NEAR( injected( q ), rate * dt * 0.31 * 0.25, 1.0e-14 );
    EXPECT_NEAR( q( 3, 4, 0 ), rate * dt, 1.0e-14 );
    EXPECT_NEAR( q( 2, 4, 0 ), rate * dt * 0.7, 1.0e-14 );
}

// Each face moves toward the inflow velocity in proportion to how much of
// its control volume the box covers
TEST( InflowSource, VelocityBlendsByOverlap )
{
    CabanaFluids::InflowSource<2> source( { 0.2, 0.2 }, { 0.3, 0.3 },
                                          { 1.0, 2.0 }, 0.0 );
    auto u = apply( source, FaceI(), 0.0, 0.5, 0.5 );
    EXPECT_NEAR( u( 3, 3, 0 ), 1.0, 1.0e-14 );
    EXPECT_NEAR( u( 2, 3, 0 ), 0.75, 1.0e-14 );
    EXPECT_NEAR( u( 6, 3, 0 ), 0.5, 1.0e-14 );

    auto v = apply( source, FaceJ(), 0.5, 0.0, 0.5 );
    EXPECT_NEAR( v( 3, 3, 0 ), 2.0, 1.0e-14 );
    EXPECT_NEAR( v( 3, 2, 0 ), 1.25, 1.0e-14 );
    EXPECT_NEAR( v( 3, 6, 0 ), 0.5, 1.0e-14 );
}
