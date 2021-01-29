/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/linear_systems/trtrs.hpp>
#include <blacspp/util/type_conversions.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  ptrtrs( blacspp::Triangle uplo, Op trans, blacspp::Diagonal diag,
    int64_t N, int64_t NRHS, 
    const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
          T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {

  auto UPLO = blacspp::detail::type_string( uplo );
  auto DIAG = blacspp::detail::type_string( diag );
  auto TRANS = char( trans );
  return wrappers::ptrtrs( UPLO.c_str(), &TRANS, DIAG.c_str(), N, NRHS,
           A, IA, JA, DESCA, B, IB, JB, DESCB );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  ptrtrs( blacspp::Triangle uplo, Op trans, blacspp::Diagonal diag,
          const BlockCyclicMatrix<T>& A, BlockCyclicMatrix<T>& B ) {

  // TODO Sanity Check
  return ptrtrs( uplo, trans, diag, B.m(), B.n(), A.data(), 1, 1, A.desc(), 
                  B.data(), 1, 1, B.desc() );

}
          
}
