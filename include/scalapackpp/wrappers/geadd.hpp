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
  pgeadd( const char* TRANS, scalapack_int M, scalapack_int N, T ALPHA,
        const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
        T BETA, 
        T* C, scalapack_int IC, scalapack_int JC, const scalapack_desc& DESCC );

}
