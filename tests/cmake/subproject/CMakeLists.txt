cmake_minimum_required( VERSION 3.17 FATAL_ERROR )
project( scalapackpp_as_subproject LANGUAGES C CXX )

# Get current Git Revision
find_package( Git REQUIRED )
execute_process( 
  COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
  OUTPUT_VARIABLE GIT_REVISION
  ERROR_QUIET
)
string( STRIP "${GIT_REVISION}" GIT_REVISION )

message( STATUS "Pulling GIT_REVISION = ${GIT_REVISION}" )

include( FetchContent )
FetchContent_Declare( scalapackpp
  GIT_REPOSITORY https://github.com/wavefunction91/scalapackpp.git
  GIT_TAG        ${GIT_REVISION}
)
FetchContent_MakeAvailable( scalapackpp )

add_executable( scalapackpp_linker_test scalapackpp_linker_test.cxx )
target_link_libraries( scalapackpp_linker_test PUBLIC scalapackpp::scalapackpp )
