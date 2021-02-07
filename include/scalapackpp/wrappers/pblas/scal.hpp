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
  pscal( int64_t N, T ALPHA, T* X, int64_t IX, int64_t JX, 
         const scalapack_desc& DESCX, int64_t INCX );
         

}
}

