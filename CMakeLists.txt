cmake_minimum_required( VERSION 3.13 FATAL_ERROR )
project( scalapackpp VERSION 0.0.0 LANGUAGES C CXX )

# Append local cmake directory to find CMAKE Modules
if( CMAKE_MODULE_PATH )
  list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else()
  set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

add_subdirectory( src )
add_subdirectory( tests )
