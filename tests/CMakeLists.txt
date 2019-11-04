include( HandleCatch2 )

add_library( ut_framework ut.cxx )
target_link_libraries( ut_framework PUBLIC scalapackpp scalapackpp::catch2 )

add_executable( test_scalapackpp scatter_gather.cxx gemm.cxx trmm.cxx trsm.cxx potrf.cxx sevp.cxx )
target_link_libraries( test_scalapackpp PUBLIC ut_framework )

