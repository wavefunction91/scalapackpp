/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/linear_systems/posv.hpp>
#include <blacspp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  pposv( blacspp::Triangle uplo, int64_t N, int64_t NRHS, 
    T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
    T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {


  auto UPLO = blacspp::detail::type_string( uplo );
  return wrappers::pposv( UPLO.c_str(), N, NRHS, A, IA, JA, DESCA,
                           B, IB, JB, DESCB );
  
}

}
