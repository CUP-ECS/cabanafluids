set(HEADERS
  BodyForce.hpp
  BoundaryConditions.hpp
  InflowSource.hpp
  Interpolation.hpp
  Mesh.hpp
  ProblemManager.hpp
  SiloWriter.hpp
  Solver.hpp
  TimeIntegrator.hpp
  VelocityCorrector.hpp
  )

set(DEPENDS_ON
  Cabana::Core
  Cabana::Grid
  SILO::silo 
  MPI::MPI_CXX
  )
 
add_library(cabanafluids INTERFACE)
set_target_properties(cabanafluids PROPERTIES INTERFACE_COMPILE_FEATURES cxx_std_17)
target_link_libraries(cabanafluids INTERFACE ${DEPENDS_ON})
target_include_directories(cabanafluids 
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS cabanafluids
  EXPORT Cabana_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
