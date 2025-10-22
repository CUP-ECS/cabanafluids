# project settings
cmake_minimum_required(VERSION 3.18)

project(CabanaFluids LANGUAGES C CXX VERSION 0.2.0)

include(GNUInstallDirs)

# Be explicit about the C++ version we want.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIIRED ON)

# find dependencies
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Get the dependent packages we need
find_package(Cabana REQUIRED COMPONENTS Cabana::Core Cabana::Grid)
if( NOT Cabana_ENABLE_MPI )
  message( FATAL_ERROR "Cabana must be compiled with MPI" )
endif()
if( NOT Cabana_ENABLE_GRID )
  message( FATAL_ERROR "Cabana must be compiled with grid support" )
endif()
find_package(SILO REQUIRED)

# Get Clang Format to use for making sure the resulting code is
# properly formatted
find_package(CLANG_FORMAT 10)

#------------------------------------------------------------------------------#
# Tests and Documentation Setup
#------------------------------------------------------------------------------#
# enable unit tests
option(CabanaFluids_ENABLE_TESTING "Build tests" ON)
if(CabanaFluids_ENABLE_TESTING)
  find_package(GTest 1.10 REQUIRED)
  # Workaround for FindGTest module in CMake older than 3.20
  if(TARGET GTest::gtest)
    set(gtest_target GTest::gtest)
  elseif(TARGET GTest::GTest)
    set(gtest_target GTest::GTest)
  else()
    message(FATAL_ERROR "bug in GTest find module workaround")
  endif()
  set(TEST_HARNESS_DIR ${CMAKE_SOURCE_DIR}/cmake/test_harness)
  include(cmake/test_harness/test_harness.cmake)
  enable_testing()
endif()

# enable doxygen
find_package(Doxygen)
if(Doxygen_FOUND)
  doxygen_add_docs(doxygen core/src grid/src)
endif()


#------------------------------------------------------------------------------#
# Libraries, Examples, and Tests
#------------------------------------------------------------------------------#
# library
add_subdirectory(src)

# examples
add_subdirectory(examples)

if(CabanaFluids_ENABLE_TESTING)
  add_subdirectory(tests)
endif()

#------------------------------------------------------------------------------#
# Clang Format
#------------------------------------------------------------------------------#
# Add a target for formatting the code using Clang
if(CLANG_FORMAT_FOUND)
  file(GLOB_RECURSE FORMAT_SOURCES src/*.cpp src/*.hpp tests/*.hpp tests/*.cpp examples/*.cpp examples/*.hpp)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXECUTABLE} -i -style=file ${FORMAT_SOURCES}
    DEPENDS ${FORMAT_SOURCES})
endif()
