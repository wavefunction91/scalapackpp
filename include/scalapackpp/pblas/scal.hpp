/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/pblas/scal.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T, typename ALPHAT>
detail::enable_if_t<
  detail::scalapack_supported_v<T> and 
  std::is_convertible_v<ALPHAT,T>
>
  pscal( int64_t N, ALPHAT ALPHA, T* X, int64_t IX, int64_t JX, 
         const scalapack_desc& DESCX, int64_t INCX ) {

  const T ALPHA_t = T(ALPHA);
  wrappers::pscal( N, ALPHA_t, X, IX, JX, DESCX, INCX );

}


}

