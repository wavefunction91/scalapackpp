/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/eigenvalue_problem/syev.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  psyev( Job jobz, Uplo uplo, int64_t N,
         T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         T* W,
         T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ ) {

  auto JOBZ = char( jobz );
  auto UPLO = char( uplo );

  int64_t LWORK = -1;
  std::vector< T > WORK( 5 );

  wrappers::psyev( &JOBZ, &UPLO, N, A, IA, JA, DESCA, 
                   W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK );

  LWORK = int64_t( WORK[0] );
  WORK.resize( LWORK );

  return wrappers::psyev( &JOBZ, &UPLO, N, A, IA, JA, DESCA, 
                          W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK );

}

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  psyevd( Job jobz, Uplo uplo, int64_t N,
         T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         T* W,
         T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ ) {

  auto JOBZ = char( jobz );
  auto UPLO = char( uplo );

  int64_t LWORK  = -1;
  int64_t LIWORK = -1;

  std::vector< T > WORK( 5 );
  std::vector< internal::scalapack_int > IWORK( 5 );

  wrappers::psyevd( &JOBZ, &UPLO, N, A, IA, JA, DESCA, 
                    W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                    IWORK.data(), LIWORK );

  LWORK  = int64_t( WORK[0]  );
  LIWORK = int64_t( IWORK[0] );
  WORK.resize( LWORK );
  IWORK.resize( LIWORK );

  return wrappers::psyevd( &JOBZ, &UPLO, N, A, IA, JA, DESCA, 
                           W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                           IWORK.data(), LIWORK );
}



template <typename T>
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  psyev( Job jobz, Uplo uplo, 
         BlockCyclicMatrix<T>& A, T* W, BlockCyclicMatrix<T>& Z ) {

  // TODO Sanity check
  return psyev( jobz, uplo, A.m(), A.data(), 1, 1, A.desc(), W, Z.data(), 1, 1,
                Z.desc() );

}

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  psyevd( Job jobz, Uplo uplo, 
          BlockCyclicMatrix<T>& A, T* W, BlockCyclicMatrix<T>& Z ) {

  // TODO Sanity check
  return psyevd( jobz, uplo, A.m(), A.data(), 1, 1, A.desc(), W, Z.data(), 1, 1,
                 Z.desc() );

}



}
