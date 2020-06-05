/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>

namespace scalapackpp::wrappers {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  paxpy( scalapack_int N, T ALPHA,
    const T* X, scalapack_int IX, scalapack_int JX, const scalapack_desc& DESCX,
          T* Y, scalapack_int IY, scalapack_int JY, const scalapack_desc& DESCY );

}
