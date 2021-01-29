/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/matrix_inverse/trtri.hpp>
#include <blacspp/util/type_conversions.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  ptrtri( blacspp::Triangle uplo, blacspp::Diagonal diag, int64_t N,
          T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA ) {

  auto UPLO = blacspp::detail::type_string( uplo );
  auto DIAG = blacspp::detail::type_string( diag );

  return wrappers::ptrtri( UPLO.c_str(), DIAG.c_str(), N, A, IA, JA, DESCA );

}


template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  ptrtri( blacspp::Triangle uplo, blacspp::Diagonal diag, BlockCyclicMatrix<T>& A ) { 

  // TODO sanity check
  return ptrtri( uplo, diag, A.m(), A.data(), 1, 1, A.desc() );

}



}



