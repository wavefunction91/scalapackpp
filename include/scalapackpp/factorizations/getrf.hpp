/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/factorizations/getrf.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgetrf( int64_t M, int64_t N, T* A, int64_t IA, int64_t JA,
          const scalapack_desc& DESCA, int64_t* IPIV ) {

  auto LOCR_A = local_row_from_desc( DESCA[internal::_M_A], DESCA );

  std::vector<internal::scalapack_int> _IPIV( LOCR_A + DESCA[internal::_MB_A] );

  auto INFO = wrappers::pgetrf( M, N, A, IA, JA, DESCA, _IPIV.data() );

  for( int64_t i = 0; i < _IPIV.size(); ++i ) IPIV[i] = _IPIV[i];

  return INFO;

}

template <typename T>
detail::enable_if_scalapack_supported_t<T,int64_t>
  pgetrf( BlockCyclicMatrix<T>& A, int64_t* IPIV ) {

  return pgetrf( A.m(), A.n(), A.data(), 1, 1, A.desc(), IPIV );

}

}
