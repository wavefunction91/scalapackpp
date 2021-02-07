/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/pblas/gemm.hpp>
#include <scalapackpp/util/type_conversions.hpp>

#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  pgemm( Op transa, Op transb,
         int64_t M, int64_t N, int64_t K, detail::type_identity_t<T> ALPHA, 
         const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         const T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB,
         detail::type_identity_t<T> BETA,
         T* C, int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {

  assert( A != C );
  assert( B != C );

  auto TRANSA = char( transa );
  auto TRANSB = char( transb );

  wrappers::pgemm( &TRANSA, &TRANSB, M, N, K, ALPHA, A, IA, JA,
                   DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC ); 

}


template <typename T>
detail::enable_if_scalapack_supported_t<T>
  pgemm( Op transa, Op transb, 
         detail::type_identity_t<T> ALPHA, const BlockCyclicMatrix<T>& A, const BlockCyclicMatrix<T>& B,
         detail::type_identity_t<T> BETA,   BlockCyclicMatrix<T>& C) {


  // TODO SANITY CHECK A/B/C/transa/b

  int64_t _M = C.m();
  int64_t _N = C.n();
  int64_t _K = transa == Op::NoTrans ? A.n() : A.m();

  pgemm( transa, transb, _M, _N, _K, ALPHA, 
         A.data(), 1, 1, A.desc(), B.data(), 1, 1, B.desc(),
         BETA, C.data(), 1, 1, C.desc() );

}

}
