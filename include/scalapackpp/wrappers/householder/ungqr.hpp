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
detail::enable_if_scalapack_complex_supported_t<T, int64_t>
  pungqr( int64_t M, int64_t N, int64_t K, T* A, int64_t IA, int64_t JA,
          const scalapack_desc& DESCA, const T* TAU, T* WORK, int64_t LWORK ); 

}
}

