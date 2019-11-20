#pragma once
#include <scalapackpp/block_cyclic.hpp>
#include <blacspp/types.hpp>

namespace scalapackpp {

template <typename T, typename ValT>
std::enable_if_t< std::is_convertible_v<ValT,T> >
  fill_triangle( const BlockCyclicDist2D& mat_dist, blacspp::Triangle uplo,
  scalapack_int M, scalapack_int N, T* A, scalapack_int LDA, const ValT& val, 
  bool include_diag = false ) {

  if( uplo == blacspp::Triangle::Upper ) {

    for( auto i = 0;                      i < M; ++i )
    for( auto j = include_diag ? i : i+1; j < N; ++j )
    if( mat_dist.i_own( i, j ) ) {
      auto [ I, J ] = mat_dist.local_indx( i, j );
      A[ I + J*LDA ]     = val;
    }

  } else {

    for( auto j = 0;                      j < N; ++j )
    for( auto i = include_diag ? j : j+1; i < M; ++i )
    if( mat_dist.i_own( i, j ) ) {
      auto [ I, J ] = mat_dist.local_indx( i, j );
      A[ I + J*LDA ]     = val;
    }

  }

}
  

}
