/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/util/sfinae.hpp>

namespace scalapackpp::wrappers {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  pgemm( const char* TRANSA, const char* TRANSB,
         scalapack_int M, scalapack_int N, scalapack_int K, T ALPHA, 
         const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         const T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
         T BETA,
         T* C, scalapack_int IC, scalapack_int JC, const scalapack_desc& DESCC );


}

