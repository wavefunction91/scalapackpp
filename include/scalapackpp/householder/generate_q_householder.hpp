/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/householder/orgqr.hpp>
#include <scalapackpp/householder/ungqr.hpp>

namespace scalapackpp {

template <
  typename T,
  detail::enable_if_scalapack_real_supported_t<T,bool> = true
>
int64_t
  generate_q_householder( int64_t M, int64_t N, int64_t K, T* A, int64_t IA, int64_t JA, 
                          const scalapack_desc& DESCA, const T* TAU ) {

  return porgqr( M, N, K, A, IA, JA, DESCA, TAU );

}

template <
  typename T,
  detail::enable_if_scalapack_complex_supported_t<T,bool> = true
>
int64_t
  generate_q_householder( int64_t M, int64_t N, int64_t K, T* A, int64_t IA, int64_t JA, 
                          const scalapack_desc& DESCA, const T* TAU ) {

  return pungqr( M, N, K, A, IA, JA, DESCA, TAU );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T,int64_t>
  generate_q_householder( int64_t K, BlockCyclicMatrix<T>& A, const T* TAU ) {

  // TODO Sanity Check
  return generate_q_householder( A.m(), A.n(), K, A.data(), 1, 1, A.desc(), TAU );

}


}
