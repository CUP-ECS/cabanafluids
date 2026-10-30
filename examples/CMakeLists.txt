add_executable( advection advection.cpp )
target_link_libraries( advection PRIVATE cabanafluids )
install( TARGETS advection DESTINATION ${CMAKE_INSTALL_BINDIR} )
