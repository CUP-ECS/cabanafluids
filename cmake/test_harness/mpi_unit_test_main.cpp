/****************************************************************************
 * Copyright (c) 2018-2023 by the Cabana authors                            *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the Cabana library. Cabana is distributed under a   *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include <gtest/gtest.h>

#include <Cabana_Grid_Config.hpp>
#include <Kokkos_Core.hpp>

#ifdef Cabana_ENABLE_HYPRE
#include <HYPRE_utilities.h>
#endif

#include <mpi.h>

int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    Kokkos::initialize( argc, argv );
#ifdef Cabana_ENABLE_HYPRE
    HYPRE_Init();
#endif
    ::testing::InitGoogleTest( &argc, argv );
    int return_val = RUN_ALL_TESTS();
#ifdef Cabana_ENABLE_HYPRE
    HYPRE_Finalize();
#endif
    Kokkos::finalize();
    MPI_Finalize();
    return return_val;
}
