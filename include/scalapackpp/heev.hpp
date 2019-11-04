#pragma once
#include <scalapackpp/wrappers/heev.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, scalapack_int>
  pheev( VectorFlag jobz, blacspp::Triangle uplo, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         detail::real_t<T>* W,
         T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ ) {

  using real_t = detail::real_t<T>;

  auto JOBZ = detail::type_string( jobz );
  auto UPLO = blacspp::detail::type_string( uplo );

  scalapack_int LWORK = -1;
  scalapack_int LRWORK = -1;
  std::vector< T > WORK( 5 );
  std::vector< real_t > RWORK( 5 );

  wrappers::pheev( JOBZ.c_str(), UPLO.c_str(), N, A, IA, JA, DESCA, 
                   W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                   RWORK.data(), LRWORK );

  LWORK  = scalapack_int( std::real(WORK[0]) );
  LRWORK = scalapack_int( RWORK[0] );
  WORK.resize( LWORK );
  RWORK.resize( 2*LRWORK );

  return wrappers::pheev( JOBZ.c_str(), UPLO.c_str(), N, A, IA, JA, DESCA, 
                          W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                          RWORK.data(), LRWORK );

}

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, scalapack_int>
  pheevd( VectorFlag jobz, blacspp::Triangle uplo, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         detail::real_t<T>* W,
         T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ ) {

  using real_t = detail::real_t<T>;

  auto JOBZ = detail::type_string( jobz );
  auto UPLO = blacspp::detail::type_string( uplo );

  scalapack_int LWORK = -1;
  scalapack_int LRWORK = -1;
  scalapack_int LIWORK = -1;
  std::vector< T > WORK( 5 );
  std::vector< real_t > RWORK( 5 );
  std::vector< scalapack_int > IWORK( 5 );

  wrappers::pheevd( JOBZ.c_str(), UPLO.c_str(), N, A, IA, JA, DESCA, 
                    W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                    RWORK.data(), LRWORK, IWORK.data(), LIWORK );

  LWORK  = scalapack_int( std::real(WORK[0]) );
  LRWORK = scalapack_int( RWORK[0] );
  LIWORK = scalapack_int( IWORK[0] );
  WORK.resize( LWORK );
  RWORK.resize( 2*LRWORK );
  IWORK.resize( LIWORK );

  return wrappers::pheevd( JOBZ.c_str(), UPLO.c_str(), N, A, IA, JA, DESCA, 
                           W, Z, IZ, JZ, DESCZ, WORK.data(), LWORK,
                           RWORK.data(), LRWORK, IWORK.data(), LIWORK );
}

}
