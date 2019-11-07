find_package( blacspp REQUIRED )
find_package( ScaLAPACK REQUIRED )

set( SCALAPACKPP_SRC
     syev.cxx 
     heev.cxx 
     gemm.cxx 
     descinit.cxx 
     gemr2d.cxx 
     type_conversions.cxx 
     trmm.cxx 
     trsm.cxx 
     information.cxx 
     potrf.cxx 
     block_cyclic.cxx
     geadd.cxx
)


add_library( scalapackpp ${SCALAPACKPP_SRC})
target_link_libraries( scalapackpp PUBLIC blacspp::blacspp ScaLAPACK::scalapack )

target_include_directories( scalapackpp PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  #$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
