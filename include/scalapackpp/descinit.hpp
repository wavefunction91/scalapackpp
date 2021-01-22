/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/descinit.hpp>
#include <blacspp/grid.hpp>

#include <tuple>

namespace scalapackpp {


std::pair< scalapack_desc, int64_t > 
inline descinit( 
  const blacspp::Grid& grid,
  int64_t M, int64_t N, int64_t MB, int64_t NB,
  int64_t ISRC, int64_t JSRC, int64_t LDD 
) {

  return wrappers::descinit( M, N, MB, NB, ISRC, JSRC, grid.context(), LDD );

}

scalapack_desc
inline descinit_noerror(
  const blacspp::Grid& grid,
  int64_t M, int64_t N, int64_t MB, int64_t NB,
  int64_t ISRC, int64_t JSRC, int64_t LDD 
) {

  scalapack_desc desc;
  int64_t        info;
  std::tie(desc,info) = descinit( grid, M, N, MB, NB, ISRC, JSRC, LDD );
  (void)info;
  return desc;

}

}
