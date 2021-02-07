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

template <
  typename T,
  detail::enable_if_scalapack_real_supported_t<T, bool> = true
>
int64_t
  pgeqpf( int64_t M, int64_t N, T* A, int64_t IA, int64_t JA,
          const scalapack_desc& DESCA, internal::scalapack_int* IPIV, 
          T* TAU, T* WORK, int64_t LWORK ); 

template <
  typename T,
  detail::enable_if_scalapack_complex_supported_t<T, bool> = true
>
int64_t
  pgeqpf( int64_t M, int64_t N, T* A, int64_t IA, int64_t JA,
          const scalapack_desc& DESCA, internal::scalapack_int* IPIV, 
          T* TAU, T* WORK, int64_t LWORK, detail::real_t<T>* RWORK,
          int64_t LRWORK ); 

}
}
