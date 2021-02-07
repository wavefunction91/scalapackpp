/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/block_cyclic.hpp>
#include <blacspp/types.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
void
  fill_triangle( const BlockCyclicDist2D& mat_dist, Uplo uplo,
  int64_t M, int64_t N, T* A, int64_t LDA, const detail::type_identity_t<T>& val, 
  bool include_diag = false ) {

  if( uplo == Uplo::Upper ) {

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
  
template <typename T>
void fill_triangle( Uplo uplo, BlockCyclicMatrix<T>& A,
                    const detail::type_identity_t<T>& val, 
                    bool include_diag = false ) {

  fill_triangle( A.dist(), uplo, A.m(), A.n(), A.data(), A.m_local(), val,
                 include_diag );

}

}
