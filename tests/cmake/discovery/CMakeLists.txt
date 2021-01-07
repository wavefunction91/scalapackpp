cmake_minimum_required( VERSION 3.17 FATAL_ERROR )
project( scalapackpp_as_subproject LANGUAGES C CXX )

find_package( scalapackpp REQUIRED )
add_executable( scalapackpp_linker_test scalapackpp_linker_test.cxx )
target_link_libraries( scalapackpp_linker_test PUBLIC scalapackpp::scalapackpp )
