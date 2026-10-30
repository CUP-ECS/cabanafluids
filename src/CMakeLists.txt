set(DEPENDS_ON
  Cabana::Core
  Cabana::Grid
  SILO::silo
  MPI::MPI_CXX
  )

add_library( cabanafluids Mesh.cpp )
target_include_directories( cabanafluids PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries( cabanafluids PUBLIC ${DEPENDS_ON} )
