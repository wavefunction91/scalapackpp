#
# a simple c++ wrapper for scalapack along with minimal extra functionality to 
# aid the the high-level development of distributed memory linear algebra.
# copyright (c) 2016-2018 david williams-young
#
# this program is free software: you can redistribute it and/or modify
# it under the terms of the gnu general public license as published by
# the free software foundation, either version 3 of the license, or
# (at your option) any later version.
#
# this program is distributed in the hope that it will be useful,
# but without any warranty; without even the implied warranty of
# merchantability or fitness for a particular purpose.  see the
# gnu general public license for more details.
#
# you should have received a copy of the gnu general public license
# along with this program.  if not, see <http://www.gnu.org/licenses/>.
#
#

# CQ Catch2 TARGET
add_library( scalapackpp::catch2 INTERFACE IMPORTED )

# Try to find Catch2
find_package( Catch2 QUIET )

if( TARGET Catch2::Catch2)
  
  message(STATUS "Found Catch2!" )
  target_link_libraries( scalapackpp::catch2 INTERFACE Catch2::Catch2 )

else()

  # Pull Catch2
  message(STATUS "Could not find Catch2! Building..." )
  include( FetchContent )
  FetchContent_Declare( catch2
    GIT_REPOSITORY      https://github.com/catchorg/Catch2.git
    GIT_TAG             v2.13.3 
    UPDATE_DISCONNECTED 1
  )

  FetchContent_MakeAvailable( catch2 )
  target_link_libraries( scalapackpp::catch2 INTERFACE Catch2 )

endif()
