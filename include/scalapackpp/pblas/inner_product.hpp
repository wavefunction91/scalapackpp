/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/pblas/dot.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T,T>
  pdot( int64_t N, const T* X, int64_t IX, int64_t JX, const scalapack_desc& DESCX, 
        int64_t INCX, const T* Y, int64_t IY, int64_t JY, const scalapack_desc& DESCY,
        int64_t INCY ) {

  T inner;
  wrappers::pdot( N, &inner, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY );
  return inner;

}

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T,T>
  pdotu( int64_t N, const T* X, int64_t IX, int64_t JX, const scalapack_desc& DESCX, 
        int64_t INCX, const T* Y, int64_t IY, int64_t JY, const scalapack_desc& DESCY,
        int64_t INCY ) {

  T inner;
  wrappers::pdotu( N, &inner, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY );
  return inner;

}

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T,T>
  pdotc( int64_t N, const T* X, int64_t IX, int64_t JX, const scalapack_desc& DESCX, 
        int64_t INCX, const T* Y, int64_t IY, int64_t JY, const scalapack_desc& DESCY,
        int64_t INCY ) {

  T inner;
  wrappers::pdotc( N, &inner, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY );
  return inner;

}

template <typename T>
detail::enable_if_scalapack_real_supported_t<T,T>
  inner_product( int64_t N, const T* X, int64_t IX, int64_t JX, 
        const scalapack_desc& DESCX, int64_t INCX, const T* Y, 
        int64_t IY, int64_t JY, const scalapack_desc& DESCY,
        int64_t INCY ) {

  return pdot( N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY );

}

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, T>
  inner_product( int64_t N, const T* X, int64_t IX, int64_t JX, 
        const scalapack_desc& DESCX, int64_t INCX, const T* Y, 
        int64_t IY, int64_t JY, const scalapack_desc& DESCY,
        int64_t INCY ) {

  return pdotc( N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY );

}


}


