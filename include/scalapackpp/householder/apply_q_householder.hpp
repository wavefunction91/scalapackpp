/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/householder/ormqr.hpp>
#include <scalapackpp/householder/unmqr.hpp>

namespace scalapackpp {

template <
  typename T,
  detail::enable_if_scalapack_real_supported_t<T,bool> = true
>
int64_t
  apply_q_householder( Side side, Op trans, int64_t M, int64_t N, 
          const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, const T* TAU,
          T* C,       int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {

  return pormqr( side, trans, M, N, A, IA, JA, DESCA, TAU, C, IC, JC, DESCC );

}

template <
  typename T,
  detail::enable_if_scalapack_complex_supported_t<T,bool> = true
>
int64_t
  apply_q_householder( Side side, Op trans, int64_t M, int64_t N, 
          const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, const T* TAU,
          T* C,       int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {

  return punmqr( side, trans, M, N, A, IA, JA, DESCA, TAU, C, IC, JC, DESCC );

}


}
