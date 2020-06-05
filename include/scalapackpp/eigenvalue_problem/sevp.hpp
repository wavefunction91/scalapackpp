/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/eigenvalue_problem/syev.hpp>
#include <scalapackpp/eigenvalue_problem/heev.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  hereig( VectorFlag jobz, blacspp::Triangle uplo, scalapack_int N,
          T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          detail::real_t<T>* W,
          T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ ) {


  if constexpr ( detail::scalapack_real_supported_v<T> )
    return psyev( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );
  else
    return pheev( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );

}


template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  hereigd( VectorFlag jobz, blacspp::Triangle uplo, scalapack_int N,
           T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
           detail::real_t<T>* W,
           T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ ) {


  if constexpr ( detail::scalapack_real_supported_v<T> )
    return psyevd( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );
  else
    return pheevd( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );

}


}
