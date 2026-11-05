/**
 * @file
 * @author Patrick Bridges <patrickb@unm.edu>
 *
 * @section DESCRIPTION
 * Inflow Sources Conditions for CabanaFluids
 */

#ifndef CABANAFLUIDS_INFLOWSOURCE_HPP
#define CABANAFLUIDS_INFLOWSOURCE_HPP

#ifndef DEBUG
#define DEBUG 0
#endif

// Include Statements
#include <Cabana_Core.hpp>
#include <Cabana_Grid.hpp>
#include <Kokkos_Core.hpp>

#include <Mesh.hpp>

namespace CabanaFluids
{
/**
 * @struct InflowSource
 * @brief Struct that applies the specified Inflow source to the underlying mesh
 */
template <std::size_t NumSpaceDim>
struct InflowSource;

template <>
struct InflowSource<2>
{
    template <class ArrayType>
    KOKKOS_INLINE_FUNCTION void operator()( Cabana::Grid::Cell, ArrayType& q,
                                            int i, int j, double x, double y,
                                            double delta_t, double v ) const
    {
        /* The inflow is a *rate*, interpolated across the volume of
         * intersecting cells, rather than a quantity forced at the location
         * as in the original incremental fluids solver */
        q( i, j, 0 ) += overlap( x, y, v ) * _rate * delta_t;
    }

    template <class ArrayType>
    KOKKOS_INLINE_FUNCTION void
    operator()( Cabana::Grid::Face<Cabana::Grid::Dim::I>, ArrayType& ux, int i,
                int j, double x, double y, [[maybe_unused]] double delta_t,
                double v ) const
    {
        double f = overlap( x, y, v );
        ux( i, j, 0 ) = ( 1.0 - f ) * ux( i, j, 0 ) + f * _velocity[0];
    }

    template <class ArrayType>
    KOKKOS_INLINE_FUNCTION void
    operator()( Cabana::Grid::Face<Cabana::Grid::Dim::J>, ArrayType& uy, int i,
                int j, double x, double y, [[maybe_unused]] double delta_t,
                double v ) const
    {
        double f = overlap( x, y, v );
        uy( i, j, 0 ) = ( 1.0 - f ) * uy( i, j, 0 ) + f * _velocity[1];
    }

    /* Fraction of the square control volume of area v centered at (x, y)
     * that the bounding box covers. Any box of positive area covers part of
     * at least one cell, however small the box. */
    KOKKOS_INLINE_FUNCTION double overlap( double x, double y, double v ) const
    {
        double half = 0.5 * sqrt( v );
        double ox = fmin( x + half, _bounding_box[2] ) -
                    fmax( x - half, _bounding_box[0] );
        double oy = fmin( y + half, _bounding_box[3] ) -
                    fmax( y - half, _bounding_box[1] );
        return ( ox > 0.0 && oy > 0.0 ) ? ox * oy / v : 0.0;
    }

    InflowSource( std::array<double, 2> location, std::array<double, 2> size,
                  std::array<double, 2> velocity, double rate )
        : _rate( rate )
    {
        _bounding_box[0] = location[0];
        _bounding_box[1] = location[1];
        _bounding_box[2] = location[0] + size[0];
        _bounding_box[3] = location[1] + size[1];
        _velocity[0] = velocity[0];
        _velocity[1] = velocity[1];
    }

    Kokkos::Array<double, 4>
        _bounding_box; /**< Bounding box on which inflow happens */
    double _rate;      /**< Quantity of advected material added per unit
                            area per unit time at the inflow location */
    Kokkos::Array<double, 2> _velocity; /**< Velocity toward which the inflow
                                           blends covered faces */
};

} // namespace CabanaFluids

#endif
