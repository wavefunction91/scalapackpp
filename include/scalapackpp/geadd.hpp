/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/geadd.hpp>
#include <scalapackpp/util/type_conversions.hpp>

#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  pgeadd( Op trans, int64_t M, int64_t N, detail::type_identity_t<T> ALPHA,
        const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
        detail::type_identity_t<T> BETA, 
        T* C, int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {

  assert( A != C );

  auto TRANS = char( trans );

  wrappers::pgeadd( &TRANS, M, N, ALPHA, A, IA, JA, DESCA, BETA,
                    C, IC, JC, DESCC );

}


template <typename T>
detail::enable_if_scalapack_supported_t<T>
  pgeadd( Op trans, detail::type_identity_t<T> ALPHA, const BlockCyclicMatrix<T>& A,
          detail::type_identity_t<T> BETA, BlockCyclicMatrix<T>& C ) {

  // TODO Sanity Check A/C/trans
  pgeadd( trans, C.m(), C.n(), ALPHA, A.data(), 1, 1, A.desc(), BETA,
          C.data(), 1, 1, C.desc() );

}

}
