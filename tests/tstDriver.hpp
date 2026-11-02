#ifndef _TSTDRIVER_HPP_
#define _TSTDRIVER_HPP_

#include "gtest/gtest.h"

#include <Kokkos_Core.hpp>

template <class E, class M>
struct DeviceType
{
    using ExecutionSpace = E;
    using MemorySpace = M;
};

// The test harness builds one executable per Kokkos backend and defines
// TEST_EXECSPACE and TEST_MEMSPACE for it
using MeshDeviceTypes =
    ::testing::Types<DeviceType<TEST_EXECSPACE, TEST_MEMSPACE>>;

#endif //_TSTDRIVER_HPP_
