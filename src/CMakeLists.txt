find_package( blacspp REQUIRED )



add_library( scalapackpp syev.cxx heev.cxx gemm.cxx descinit.cxx gemr2d.cxx type_conversions.cxx trmm.cxx trsm.cxx information.cxx potrf.cxx block_cyclic.cxx )
target_link_libraries( scalapackpp PUBLIC blacspp::blacspp )

target_include_directories( scalapackpp PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  #$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
