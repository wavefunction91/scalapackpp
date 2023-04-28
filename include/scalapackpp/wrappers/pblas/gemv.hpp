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
  pgemv( const char* TRANS,
         int64_t M, int64_t N, T ALPHA, 
         const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         const T* X, int64_t IX, int64_t JX, const scalapack_desc& DESCX, int64_t INCX,
         T BETA,
         T* Y, int64_t IY, int64_t JY, const scalapack_desc& DESCY, int64_t INCY );


}
}

