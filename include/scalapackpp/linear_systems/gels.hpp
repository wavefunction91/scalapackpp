/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/linear_systems/gels.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgels( TransposeFlag trans, int64_t M, int64_t N, int64_t NRHS, 
    T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
    T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {

  auto TRANS = detail::type_string( trans );

  int64_t LWORK = -1;
  std::vector< T > WORK(5);

  wrappers::pgels( TRANS.c_str(), M, N, NRHS, A, IA, JA, DESCA, 
                   B, IB, JB, DESCB, WORK.data(), LWORK );

  LWORK = int64_t( std::real(WORK[0]) );
  WORK.resize( LWORK );

  return wrappers::pgels( TRANS.c_str(), M, N, NRHS, A, IA, JA, DESCA, 
                          B, IB, JB, DESCB, WORK.data(), LWORK );

}

}
