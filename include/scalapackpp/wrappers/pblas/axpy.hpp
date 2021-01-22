/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>

namespace scalapackpp {
namespace wrappers    {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  paxpy( int64_t N, T ALPHA,
    const T* X, int64_t IX, int64_t JX, const scalapack_desc& DESCX,
          T* Y, int64_t IY, int64_t JY, const scalapack_desc& DESCY );

}
}
