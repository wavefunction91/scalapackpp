/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/factorizations/getrf.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  pgetrf( scalapack_int M, scalapack_int N, T* A, scalapack_int IA, scalapack_int JA,
          const scalapack_desc& DESCA, scalapack_int* IPIV ) {

  return wrappers::pgetrf( M, N, A, IA, JA, DESCA, IPIV );

}

}
