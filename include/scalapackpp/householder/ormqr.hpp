/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/householder/ormqr.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  pormqr( SideFlag side, TransposeFlag trans, int64_t M, int64_t N, 
          const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, const T* TAU,
          T* C,       int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {


  auto SIDE  = detail::type_string( side  );
  auto TRANS = detail::type_srting( trans );

  int64_t LWORK = -1;
  std::vector<T> WORK(5);

  wrappers::pormqr( SIDE.c_str(), TRANS.c_str(), M, N, A, IA, JA, DESCA, TAU,
                    C, IC, JC, DESCC, WORK.data(), LWORK );

  LWORK = int64_t( WORK[0] );
  WORK.resize( LWORK );

  return wrappers::pormqr( SIDE.c_str(), TRANS.c_str(), M, N, A, IA, JA, DESCA, TAU,
                           C, IC, JC, DESCC, WORK.data(), LWORK );

}

}
