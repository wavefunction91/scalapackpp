/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/linear_systems/getrs.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/wrappers/support.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgetrs( TransposeFlag trans, int64_t N, int64_t NRHS, 
    const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
    int64_t* IPIV,
    T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {

  auto TRANS = detail::type_string( trans );

  const auto ICXT_A = DESCA[internal::_CTXT_A];
  const auto M_A    = DESCA[internal::_M_A];
  const auto MB_A   = DESCA[internal::_MB_A];
  const auto RSRC_A = DESCA[internal::_RSRC_A];

  auto grid_dim = blacspp::wrappers::grid_info( ICXT_A );
  const auto LOCR_A = numroc( M_A, MB_A, RSRC_A, grid_dim.my_row, grid_dim.np_row );

  std::vector<internal::scalapack_int> _IPIV( LOCR_A + MB_A );
  for( int64_t i = 0; i < _IPIV.size(); ++i )
    _IPIV[i] = IPIV[i] ;

  return wrappers::pgetrs( TRANS.c_str(), N, NRHS, A, IA, JA, DESCA, 
    _IPIV.data(), B, IB, JB, DESCB );

}
          
}
