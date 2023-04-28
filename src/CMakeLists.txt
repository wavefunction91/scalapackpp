if( NOT TARGET blacspp::blacspp )

  find_package( blacspp QUIET ) # blacspp links to ScaLAPACK
  if( NOT blacspp_FOUND )

    set( BLACSPP_ENABLE_ILP64 ${SCALAPACKPP_ENABLE_ILP64} CACHE BOOL "" )
    set( BLACSPP_FORCE_ILP64  ${SCALAPACKPP_FORCE_ILP64}  CACHE BOOL "" )
    include( FetchContent )
    FetchContent_Declare( blacspp
      GIT_REPOSITORY https://github.com/wavefunction91/blacspp.git
      GIT_TAG        f5cd6150d9b18790a17a5de667c44e300eafb3ad
    )
    FetchContent_MakeAvailable( blacspp )

  else()

    if( SCALAPACKPP_FORCE_ILP64 AND ScaLAPACK_IS_LP64 )
      message( FATAL_ERROR 
               "ScaLAPACK++ Forced ILP64 - ScaLAPACK in BLACS++ is LP64" )
    endif()

  endif()

endif()

set( ScaLAPACKPP_SRC
     syev.cxx 
     heev.cxx 
     sygst.cxx 
     hegst.cxx 
     gemm.cxx 
     gemv.cxx 
     symm.cxx 
     symv.cxx 
     hemv.cxx 
     hemm.cxx 
     descinit.cxx 
     gemr2d.cxx 
     gesv.cxx
     gels.cxx
     getrf.cxx
     geqrf.cxx
     geqpf.cxx
     getrs.cxx
     getri.cxx
     orgqr.cxx
     ungqr.cxx
     ormqr.cxx
     unmqr.cxx
     gesvd.cxx
     type_conversions.cxx 
     trmm.cxx 
     trsm.cxx 
     trtri.cxx
     trtrs.cxx
     information.cxx 
     posv.cxx
     potrf.cxx 
     potrs.cxx
     potri.cxx
     block_cyclic.cxx
     geadd.cxx
     tradd.cxx
     lange.cxx
     lansy.cxx
     lanhe.cxx
     scal.cxx
     axpy.cxx
     copy.cxx
     dot.cxx
     nrm2.cxx
     lacpy.cxx
     lascl.cxx
)


add_library( scalapackpp ${ScaLAPACKPP_SRC})
target_link_libraries( scalapackpp PUBLIC blacspp::blacspp )

target_include_directories( scalapackpp PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  #$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
)


# INSTALL rules
add_library( scalapackpp::scalapackpp ALIAS scalapackpp )

include( GNUInstallDirs )

# TARGETS

# Check if blacspp is imported (e.g. not a subproject)
get_target_property( BLACSPP_IS_IMPORTED blacspp::blacspp IMPORTED )

# Export blacspp if not imported
set( BLASPP_EXPORT_TARGET )
if( NOT BLACSPP_IS_IMPORTED )
  message( STATUS "SCALAPACKPP WILL EXPORT BLACSPP" )
  set( BLACSPP_EXPORT_TARGET blacspp )
endif()

install( TARGETS scalapackpp ${BLACSPP_EXPORT_TARGET}
  EXPORT scalapackpp-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set_target_properties( scalapackpp PROPERTIES EXPORT_NAME scalapackpp )

# Export build tree
export(EXPORT scalapackpp-targets
      NAMESPACE scalapackpp::
      FILE "${PROJECT_BINARY_DIR}/scalapackpp-targets.cmake")

# Install Headers
install( 
  DIRECTORY   ${PROJECT_SOURCE_DIR}/include
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  FILES_MATCHING PATTERN "*.hpp"
)


# Export target to scripe
set( INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/scalapackpp )
install( EXPORT scalapackpp-targets
  FILE         scalapackpp-targets.cmake
  NAMESPACE    scalapackpp::
  DESTINATION  ${INSTALL_CONFIGDIR}
)


# Create config-version.cmake file
include( CMakePackageConfigHelpers )
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/scalapackpp-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)



# Setup scalapackpp-config.cmake
configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/scalapackpp-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/scalapackpp-config.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Install misc cmake modules
install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/scalapackpp-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/scalapackpp-config-version.cmake
  DESTINATION ${INSTALL_CONFIGDIR}
)


