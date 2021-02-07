/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/householder/orgqr.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  porgqr( int64_t M, int64_t N, int64_t K, T* A, int64_t IA, int64_t JA, 
          const scalapack_desc& DESCA, const T* TAU ) {

  int64_t LWORK = -1;
  std::vector<T> WORK(5);

  wrappers::porgqr( M, N, K, A, IA, JA, DESCA, TAU, WORK.data(), LWORK );

  LWORK = int64_t( WORK[0] );
  WORK.resize( LWORK );

  return wrappers::porgqr( M, N, K, A, IA, JA, DESCA, TAU, WORK.data(), LWORK );

}

}
