/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/matrix_inverse/trtri.hpp>
#include <blacspp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, scalapack_int>
  ptrtri( blacspp::Triangle uplo, blacspp::Diagonal diag, scalapack_int N,
          T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA ) {

  auto UPLO = blacspp::detail::type_string( uplo );
  auto DIAG = blacspp::detail::type_string( diag );

  return wrappers::ptrtri( UPLO.c_str(), DIAG.c_str(), N, A, IA, JA, DESCA );

}

/*
template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, scalapack_int>
  ptrtri( blacspp::Triangle uplo, blacspp::Diagonal diag, scalapack_int N,
          T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA ) {

  static_assert(false, "COMPLEX PXTRTRI SEEMS TO BE BROKEN... FILE A BUG REPORT WITH THE DEVELOPERS");

  auto UPLO = blacspp::detail::type_string( uplo );
  auto DIAG = blacspp::detail::type_string( diag );

  return wrappers::ptrtri( UPLO.c_str(), DIAG.c_str(), N, A, IA, JA, DESCA );

}
*/


}



