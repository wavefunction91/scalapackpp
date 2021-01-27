/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/linear_systems/gesv.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgesv( int64_t N, int64_t NRHS, 
    T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
    int64_t* IPIV,
    T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {

  auto LOCR_A = local_row_from_desc( DESCA[internal::_M_A], DESCA );

  std::vector<internal::scalapack_int> _IPIV( LOCR_A + DESCA[internal::_MB_A] );

  auto INFO = wrappers::pgesv( N, NRHS, A, IA, JA, DESCA, _IPIV.data(), 
                               B, IB, JB, DESCB );

  for( int64_t i = 0; i < _IPIV.size(); ++i ) IPIV[i] = _IPIV[i];

  return INFO;
}
          
}
