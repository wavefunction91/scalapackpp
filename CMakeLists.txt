cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
project( scalapackpp VERSION 0.0.0 LANGUAGES C CXX )

# Append local cmake directory to find CMAKE Modules
if( CMAKE_MODULE_PATH )
  list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else()
  set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

add_subdirectory( src )

if(NOT DEFINED SCALAPACKPP_ENABLE_TESTS )
  set( SCALAPACKPP_ENABLE_TESTS ON )
endif()

if( SCALAPACKPP_ENABLE_TESTS )
  add_subdirectory( tests )
endif()
