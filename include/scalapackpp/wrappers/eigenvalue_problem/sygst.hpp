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
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  psygst( int64_t IBTYPE, const char* UPLO, int64_t N,
                T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
          const T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB,
                T* SCALE );

}
}
