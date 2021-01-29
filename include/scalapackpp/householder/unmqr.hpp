/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/householder/unmqr.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, int64_t>
  punmqr( Side side, Op trans, int64_t M, int64_t N, 
          const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, const T* TAU,
          T* C,       int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {


  auto SIDE  = char( side  );
  auto TRANS = char( trans );

  int64_t LWORK = -1;
  std::vector<T> WORK(5);

  wrappers::punmqr( &SIDE, &TRANS, M, N, A, IA, JA, DESCA, TAU,
                    C, IC, JC, DESCC, WORK.data(), LWORK );

  LWORK = int64_t( std::real(WORK[0]) );
  WORK.resize( LWORK );

  return wrappers::punmqr( &SIDE, &TRANS, M, N, A, IA, JA, DESCA, TAU,
                           C, IC, JC, DESCC, WORK.data(), LWORK );

}

}
