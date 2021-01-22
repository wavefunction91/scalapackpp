/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/geadd.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  transpose( int64_t M, int64_t N, 
    const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
          T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB 
  ) {

  pgeadd( TransposeFlag::Transpose, M, N, 1., A, IA, JA, DESCA, 
          0., B, IB, JB, DESCB );

}

template <typename T>
detail::enable_if_scalapack_real_supported_t<T>
  conj_transpose( int64_t M, int64_t N, 
    const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
          T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB 
  ) {

  transpose( M, N, A, IA, JA, DESCA, B, IB, JB, DESCB );

}
template <typename T>
detail::enable_if_scalapack_complex_supported_t<T>
  conj_transpose( int64_t M, int64_t N, 
    const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
          T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB 
  ) {

  pgeadd( TransposeFlag::ConjTranspose, M, N, 1., A, IA, JA, DESCA, 
          0., B, IB, JB, DESCB );

}
          

} 
