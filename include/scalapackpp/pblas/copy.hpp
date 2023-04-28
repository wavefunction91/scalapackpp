/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/pblas/copy.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  pcopy( int64_t N, const T* X, int64_t IX, int64_t JX, 
         const scalapack_desc& DESCX, int64_t INCX, T* Y, int64_t IY,
         int64_t JY, const scalapack_desc& DESCY, int64_t INCY ) {

  wrappers::pcopy( N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY );

}


}

