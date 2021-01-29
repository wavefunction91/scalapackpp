/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/pblas/trsm.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  ptrsm( SideFlag side, blacspp::Triangle uplo, TransposeFlag trans, 
         blacspp::Diagonal diag,
         int64_t M, int64_t N, detail::type_identity_t<T> ALPHA, 
         const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {

  auto SIDE = detail::type_string( side );
  auto UPLO = blacspp::detail::type_string( uplo );
  auto TRANS = detail::type_string( trans );
  auto DIAG  = blacspp::detail::type_string( diag );


  wrappers::ptrsm( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
                   M, N, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  ptrsm( SideFlag side, blacspp::Triangle uplo, TransposeFlag trans, 
         blacspp::Diagonal diag, detail::type_identity_t<T> ALPHA, 
         const BlockCyclicMatrix<T>& A, BlockCyclicMatrix<T>& B ) {

  // TODO Sanity check
  ptrsm( side, uplo, trans, diag, B.m(), B.n(), ALPHA, A.data(), 1, 1, A.desc(),
         B.data(), 1, 1, B.desc() );
}

}
