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
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  hereig( Job jobz, Uplo uplo, int64_t N,
          T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
          detail::real_t<T>* W,
          T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ ) {


  return psyev( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );

}

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, int64_t>
  hereig( Job jobz, Uplo uplo, int64_t N,
          T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
          detail::real_t<T>* W,
          T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ ) {


  return pheev( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );

}


template <typename T>
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  hereigd( Job jobz, Uplo uplo, int64_t N,
           T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
           detail::real_t<T>* W,
           T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ ) {


  return psyevd( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );

}

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, int64_t>
  hereigd( Job jobz, Uplo uplo, int64_t N,
           T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
           detail::real_t<T>* W,
           T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ ) {


  return pheevd( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );

}



template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  hereig( Job jobz, Uplo uplo, 
          BlockCyclicMatrix<T>& A, detail::real_t<T>* W, BlockCyclicMatrix<T>& Z ) {

  // TODO Sanity check
  return hereig( jobz, uplo, A.m(), A.data(), 1, 1, A.desc(), W, Z.data(), 1, 1,
                 Z.desc() );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  hereigd( Job jobz, Uplo uplo, 
           BlockCyclicMatrix<T>& A, detail::real_t<T>* W, BlockCyclicMatrix<T>& Z ) {

  // TODO Sanity check
  return hereigd( jobz, uplo, A.m(), A.data(), 1, 1, A.desc(), W, Z.data(), 1, 1,
                  Z.desc() );

}


}
