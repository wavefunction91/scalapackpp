/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/matrix_inverse/potri.hpp>
#include <blacspp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  ppotri( blacspp::Triangle uplo, int64_t N, T* A, int64_t IA, 
          int64_t JA, const scalapack_desc& DESCA ) {

  auto UPLO = blacspp::detail::type_string( uplo );
  return wrappers::ppotri( UPLO.c_str(), N, A, IA, JA, DESCA );

}

}
