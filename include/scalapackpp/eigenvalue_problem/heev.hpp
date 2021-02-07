/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/eigenvalue_problem/heev.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, int64_t>
  pheev( Job jobz, Uplo uplo, int64_t N,
         T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         detail::real_t<T>* W,
         T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ ) {

  using real_t = detail::real_t<T>;

  auto JOBZ = char( jobz );
  auto UPLO = char( uplo );

  int64_t LWORK = -1;
  int64_t LRWORK = -1;
  std::vector< T > WORK( 5 );
  std::vector< real_t > RWORK( 5 );

  wrappers::pheev( &JOBZ, &UPLO, N, A, IA, JA, DESCA, 
                   W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                   RWORK.data(), LRWORK );

  LWORK  = int64_t( std::real(WORK[0]) );
  LRWORK = int64_t( RWORK[0] );
  WORK.resize( LWORK );
  RWORK.resize( 2*LRWORK );

  return wrappers::pheev( &JOBZ, &UPLO, N, A, IA, JA, DESCA, 
                          W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                          RWORK.data(), LRWORK );

}

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, int64_t>
  pheevd( Job jobz, Uplo uplo, int64_t N,
         T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         detail::real_t<T>* W,
         T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ ) {

  using real_t = detail::real_t<T>;

  auto JOBZ = char( jobz );
  auto UPLO = char( uplo );

  int64_t LWORK = -1;
  int64_t LRWORK = -1;
  int64_t LIWORK = -1;
  std::vector< T > WORK( 5 );
  std::vector< real_t > RWORK( 5 );
  std::vector< internal::scalapack_int > IWORK( 5 );

  wrappers::pheevd( &JOBZ, &UPLO, N, A, IA, JA, DESCA, 
                    W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                    RWORK.data(), LRWORK, IWORK.data(), LIWORK );

  LWORK  = int64_t( std::real(WORK[0]) );
  LRWORK = int64_t( RWORK[0] );
  LIWORK = int64_t( IWORK[0] );
  WORK.resize( LWORK );
  RWORK.resize( 2*LRWORK );
  IWORK.resize( LIWORK );

  return wrappers::pheevd( &JOBZ, &UPLO, N, A, IA, JA, DESCA, 
                           W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                           RWORK.data(), LRWORK, IWORK.data(), LIWORK );
}



template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, int64_t>
  pheev( Job jobz, Uplo uplo, 
         BlockCyclicMatrix<T>& A, detail::real_t<T>* W, BlockCyclicMatrix<T>& Z ) {

  // TODO Sanity check
  return pheev( jobz, uplo, A.m(), A.data(), 1, 1, A.desc(), W, Z.data(), 1, 1,
                Z.desc() );

}

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, int64_t>
  pheevd( Job jobz, Uplo uplo, 
          BlockCyclicMatrix<T>& A, detail::real_t<T>* W, BlockCyclicMatrix<T>& Z ) {

  // TODO Sanity check
  return pheevd( jobz, uplo, A.m(), A.data(), 1, 1, A.desc(), W, Z.data(), 1, 1,
                 Z.desc() );

}

}
