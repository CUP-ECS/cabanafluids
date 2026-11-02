foreach( _test Mesh ProblemManager BoundaryConditions TimeIntegrator VelocityCorrector InflowSource )
  add_executable( tst${_test} tst${_test}.cpp )
  target_link_libraries( tst${_test} PRIVATE cabanafluids GTest::gtest )
  foreach( _np 1 2 4 )
    add_test( NAME ${_test}Tests_${_np}
              COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_np}
                      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:tst${_test}>
                      ${MPIEXEC_POSTFLAGS} )
    set_tests_properties( ${_test}Tests_${_np} PROPERTIES PROCESSORS ${_np} )
  endforeach()
endforeach()
