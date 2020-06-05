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

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, scalapack_int>
  psyev( VectorFlag jobz, blacspp::Triangle uplo, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         T* W,
         T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ ) {

  auto JOBZ = detail::type_string( jobz );
  auto UPLO = blacspp::detail::type_string( uplo );

  scalapack_int LWORK = -1;
  std::vector< T > WORK( 5 );

  wrappers::psyev( JOBZ.c_str(), UPLO.c_str(), N, A, IA, JA, DESCA, 
                   W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK );

  LWORK = scalapack_int( WORK[0] );
  WORK.resize( LWORK );

  return wrappers::psyev( JOBZ.c_str(), UPLO.c_str(), N, A, IA, JA, DESCA, 
                          W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK );

}

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, scalapack_int>
  psyevd( VectorFlag jobz, blacspp::Triangle uplo, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         T* W,
         T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ ) {

  auto JOBZ = detail::type_string( jobz );
  auto UPLO = blacspp::detail::type_string( uplo );

  scalapack_int LWORK  = -1;
  scalapack_int LIWORK = -1;

  std::vector< T > WORK( 5 );
  std::vector< scalapack_int > IWORK( 5 );

  wrappers::psyevd( JOBZ.c_str(), UPLO.c_str(), N, A, IA, JA, DESCA, 
                    W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                    IWORK.data(), LIWORK );

  LWORK  = scalapack_int( WORK[0]  );
  LIWORK = scalapack_int( IWORK[0] );
  WORK.resize( LWORK );
  IWORK.resize( LIWORK );

  return wrappers::psyevd( JOBZ.c_str(), UPLO.c_str(), N, A, IA, JA, DESCA, 
                           W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                           IWORK.data(), LIWORK );
}

}
