/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/factorizations/geqrf.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgeqrf( int64_t M, int64_t N, T* A, int64_t IA, int64_t JA,
          const scalapack_desc& DESCA, T* TAU ) {

  int64_t LWORK = -1;
  std::vector< T > WORK(5);

  wrappers::pgeqrf( M, N, A, IA, JA, DESCA, TAU, WORK.data(), LWORK );

  LWORK = int64_t( std::real(WORK[0]) );
  WORK.resize( LWORK );

  return wrappers::pgeqrf( M, N, A, IA, JA, DESCA, TAU, WORK.data(), LWORK );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T,int64_t>
  pgeqrf( BlockCyclicMatrix<T>& A, T* TAU ) {

  return pgeqrf( A.m(), A.n(), A.data(), 1, 1, A.desc(), TAU );

}

}

