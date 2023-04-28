/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/pblas/axpy.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T, typename ALPHAT>
detail::enable_if_t<
  detail::scalapack_supported_v<T> and 
  std::is_convertible_v<ALPHAT,T>
>
  paxpy( int64_t N, ALPHAT ALPHA, const T* X, int64_t IX, int64_t JX, 
         const scalapack_desc& DESCX, int64_t INCX, T* Y, int64_t IY,
         int64_t JY, const scalapack_desc& DESCY, int64_t INCY ) {

  const T ALPHA_t = T(ALPHA);
  wrappers::paxpy( N, ALPHA_t, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY );

}


}

