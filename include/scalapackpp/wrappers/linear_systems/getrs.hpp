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
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgetrs( const char* TRANS, int64_t N, int64_t NRHS, 
    const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
    const internal::scalapack_int* IPIV,
    T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB );
          
}
}
