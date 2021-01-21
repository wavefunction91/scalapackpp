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
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgetrf( int64_t M, int64_t N, T* A, int64_t IA, int64_t JA,
          const scalapack_desc& DESCA, internal::scalapack_int* IPIV );

}
