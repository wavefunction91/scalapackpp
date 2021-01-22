/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/util/sfinae.hpp>

namespace scalapackpp {
namespace wrappers    {

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, detail::real_t<T>>
  planhe( const char* NORM, const char* UPLO, int64_t N, const T* A, 
          int64_t IA, int64_t JA, const scalapack_desc& DESCA,
          detail::real_t<T>* WORK );

}
}

