if( NOT TARGET blacspp::blacspp )
  find_package( blacspp REQUIRED )
endif()
find_package( ScaLAPACK REQUIRED )

set( ScaLAPACKPP_SRC
     syev.cxx 
     heev.cxx 
     gemm.cxx 
     descinit.cxx 
     gemr2d.cxx 
     gesv.cxx
     getrf.cxx
     getrs.cxx
     getri.cxx
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
     scal.cxx
)


add_library( scalapackpp ${ScaLAPACKPP_SRC})
target_link_libraries( scalapackpp PUBLIC blacspp::blacspp ScaLAPACK::scalapack )

target_include_directories( scalapackpp PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  #$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
)


# INSTALL rules
add_library( scalapackpp::scalapackpp ALIAS scalapackpp )

include( GNUInstallDirs )

# TARGETS
install( TARGETS scalapackpp
  EXPORT scalapackpp-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set_target_properties( scalapackpp PROPERTIES EXPORT_NAME scalapackpp )

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
  ${PROJECT_SOURCE_DIR}/cmake/modules/CommonFunctions.cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules/FindBLACS.cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules/FindBLAS.cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules/FindLAPACK.cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules/FindScaLAPACK.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/scalapackpp-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/scalapackpp-config-version.cmake
  DESTINATION ${INSTALL_CONFIGDIR}
)
install( FILES  
  ${PROJECT_SOURCE_DIR}/cmake/modules/linalg-modules/FindBLIS.cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules/linalg-modules/FindIBMESSL.cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules/linalg-modules/FindIntelMKL.cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules/linalg-modules/FindOpenBLAS.cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules/linalg-modules/FindReferenceBLAS.cmake
  ${PROJECT_SOURCE_DIR}/cmake/modules/linalg-modules/FindReferenceLAPACK.cmake
  DESTINATION ${INSTALL_CONFIGDIR}/linalg-modules
)


