foreach( _test Mesh ProblemManager BoundaryConditions )
  add_executable( tst${_test} tst${_test}.cpp )
  target_link_libraries( tst${_test} PRIVATE cabanafluids GTest::gtest )
  add_test( NAME ${_test}Tests COMMAND tst${_test} )
endforeach()

# tstInflowSource does not exist yet
#add_executable( tstInflowSource tstInflowSource.cpp )
#target_link_libraries( tstInflowSource PRIVATE cabanafluids GTest::gtest )
#add_test( NAME InflowSourceTests COMMAND tstInflowSource )
