/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/block_cyclic.hpp>
#include <blacspp/types.hpp>

namespace scalapackpp {

template <typename T, typename ValT>
std::enable_if_t< std::is_convertible_v<ValT,T> >
  fill_triangle( const BlockCyclicDist2D& mat_dist, blacspp::Triangle uplo,
  int64_t M, int64_t N, T* A, int64_t LDA, const ValT& val, 
  bool include_diag = false ) {

  if( uplo == blacspp::Triangle::Upper ) {

    for( int64_t i = 0;                      i < M; ++i )
    for( int64_t j = include_diag ? i : i+1; j < N; ++j )
    if( mat_dist.i_own( i, j ) ) {
      auto [ I, J ] = mat_dist.local_indx( i, j );
      A[ I + J*LDA ]     = val;
    }

  } else {

    for( int64_t j = 0;                      j < N; ++j )
    for( int64_t i = include_diag ? j : j+1; i < M; ++i )
    if( mat_dist.i_own( i, j ) ) {
      auto [ I, J ] = mat_dist.local_indx( i, j );
      A[ I + J*LDA ]     = val;
    }

  }

}
  

}
