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


  int64_t local_row_from_desc( int64_t M, const scalapack_desc& desc );
  int64_t local_col_from_desc( int64_t N, const scalapack_desc& desc );


}
