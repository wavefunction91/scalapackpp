/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/util/type_traits.hpp>

namespace scalapackpp {
namespace wrappers    {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  pgemm( const char* TRANSA, const char* TRANSB,
         int64_t M, int64_t N, int64_t K, T ALPHA, 
         const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         const T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB,
         T BETA,
         T* C, int64_t IC, int64_t JC, const scalapack_desc& DESCC );


}
}

