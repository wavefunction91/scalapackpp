/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/lascl.hpp>
#include <scalapackpp/util/type_conversions.hpp>

#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  plascl( MatrixType type, detail::real_t<T> CFROM, 
          detail::real_t<T> CTO, int64_t M, int64_t N, 
          T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA ) {

  auto TYPE = char( type );
  wrappers::plascl( &TYPE, CFROM, CTO, M, N, A, IA, JA, DESCA );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  plascl( MatrixType type, detail::real_t<T> ALPHA, int64_t M, int64_t N, 
          T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA ) {

  plascl( type, 1., ALPHA,  M, N, A, IA, JA, DESCA );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  plascl( MatrixType type, detail::real_t<T> CFROM, 
    detail::real_t<T> CTO, BlockCyclicMatrix<T>& A ) {

  plascl( type, CFROM, CTO, A.m(), A.n(), A.data(), 1, 1, A.desc() );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  plascl( MatrixType type, detail::real_t<T> ALPHA, BlockCyclicMatrix<T>& A ) {

  plascl( type, 1., ALPHA, A );

}

}
