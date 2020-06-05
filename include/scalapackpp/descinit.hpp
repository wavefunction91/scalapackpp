/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/descinit.hpp>
#include <blacspp/grid.hpp>

namespace scalapackpp {


std::pair< scalapack_desc, scalapack_int > 
inline descinit( 
  const blacspp::Grid& grid,
  scalapack_int M, scalapack_int N, scalapack_int MB, scalapack_int NB,
  scalapack_int ISRC, scalapack_int JSRC, scalapack_int LDD 
) {

  return wrappers::descinit( M, N, MB, NB, ISRC, JSRC, grid.context(), LDD );

}

scalapack_desc
inline descinit_noerror(
  const blacspp::Grid& grid,
  scalapack_int M, scalapack_int N, scalapack_int MB, scalapack_int NB,
  scalapack_int ISRC, scalapack_int JSRC, scalapack_int LDD 
) {

  auto [ desc, info ] = descinit( grid, M, N, MB, NB, ISRC, JSRC, LDD );
  return desc;

}

}
