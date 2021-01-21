/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <blacspp/grid.hpp>

namespace scalapackpp {

  int64_t numroc( int64_t N, int64_t NB,
    int64_t IPROC, int64_t ISRC, 
    int64_t NPROC );

  std::pair<int64_t, int64_t>
    get_local_dims( const blacspp::Grid& grid,
      int64_t M, int64_t N,
      int64_t MB, int64_t NB,
      int64_t ISRC, int64_t JSRC );

/*
  std::pair< block_cyclic_coordinate, block_cyclic_coordinate >
    local_index( const blacspp::Grid& grid,
      int64_t MB, int64_t NB,
      int64_t I, int64_t J
    );

  std::pair< int64_t, int64_t >
    global_index( const blacspp::Grid& grid,
      int64_t MB, int64_t NB,
      int64_t I, int64_t J
    );
*/
}
