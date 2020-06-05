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
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  pgetrs( const char* TRANS, scalapack_int N, scalapack_int NRHS, 
    const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    const scalapack_int* IPIV,
    T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB );
          
}
