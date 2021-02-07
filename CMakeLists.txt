#
# This file is a part of scalapackspp (see LICENSE)
#
# Copyright (c) 2019-2020 David Williams-Young
# All rights reserved
#
cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
project( scalapackpp VERSION 0.0.0 LANGUAGES C CXX )

# Append local cmake directory to find CMAKE Modules
if( CMAKE_MODULE_PATH )
  list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else()
  set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()


# Configure Options
include( CMakeDependentOption )
option( SCALAPACKPP_ENABLE_ILP64 "Enable search for ILP64 ScaLAPACK bindings" OFF )
cmake_dependent_option( SCALAPACKPP_FORCE_ILP64 "Force ILP64 - Fail if not found" OFF
                        "SCALAPACKPP_ENABLE_ILP64" OFF )

add_subdirectory( src )

if(NOT DEFINED SCALAPACKPP_ENABLE_TESTS )
  set( SCALAPACKPP_ENABLE_TESTS ON )
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

if( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND SCALAPACKPP_ENABLE_TESTS AND BUILD_TESTING )
  enable_testing()
  add_subdirectory( tests )
endif()
