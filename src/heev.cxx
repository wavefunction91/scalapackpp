/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/eigenvalue_problem/heev.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pcheev_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, float* W,
         scomplex* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, scomplex* WORK, const scalapack_int *LWORK,
         float* RWORK, const scalapack_int *LRWORK, scalapack_int* INFO );

void pzheev_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, double* W,
         dcomplex* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, dcomplex* WORK, const scalapack_int *LWORK,
         double* RWORK, const scalapack_int *LRWORK, scalapack_int* INFO );

void pcheevd_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, float* W,
         scomplex* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, scomplex* WORK, const scalapack_int *LWORK,
         float* RWORK, const scalapack_int* LRWORK,
         scalapack_int* IWORK, const scalapack_int* LIWORK, scalapack_int* INFO );

void pzheevd_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, double* W,
         dcomplex* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, dcomplex* WORK, const scalapack_int *LWORK,
         double* RWORK, const scalapack_int* LRWORK,
         scalapack_int* IWORK, const scalapack_int* LIWORK, scalapack_int* INFO );

}


namespace scalapackpp::wrappers {

#define pheev_impl(F, FNAME)                                        \
template <>                                                         \
int64_t                                                             \
  pheev( const char* JOBZ, const char* UPLO, int64_t N,             \
         F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, \
         detail::real_t<F>* W,                                      \
         F* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ, \
         F* WORK, int64_t LWORK,                                    \
         detail::real_t<F>* RWORK, int64_t LRWORK ) {               \
                                                                    \
  auto _N      = detail::to_scalapack_int( N  );                    \
  auto _IA     = detail::to_scalapack_int( IA );                    \
  auto _JA     = detail::to_scalapack_int( JA );                    \
  auto _IZ     = detail::to_scalapack_int( IZ );                    \
  auto _JZ     = detail::to_scalapack_int( JZ );                    \
  auto _LWORK  = detail::to_scalapack_int( LWORK );                 \
  auto _LRWORK = detail::to_scalapack_int( LRWORK );                \
                                                                    \
  auto _DESCA = detail::to_scalapack_int( DESCA );                  \
  auto _DESCZ = detail::to_scalapack_int( DESCZ );                  \
                                                                    \
  scalapack_int INFO;                                               \
  FNAME( JOBZ, UPLO, &_N, A, &_IA, &_JA, _DESCA.data(),             \
         W, Z, &_IZ, &_JZ, _DESCZ.data(), WORK, &_LWORK,            \
         RWORK, &_LRWORK, &INFO );                                  \
  return INFO;                                                      \
                                                                    \
}

pheev_impl( scomplex, pcheev_ );
pheev_impl( dcomplex, pzheev_ );

#define pheevd_impl(F, FNAME)                                       \
template <>                                                         \
int64_t                                                             \
  pheevd( const char* JOBZ, const char* UPLO, int64_t N,            \
         F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, \
         detail::real_t<F>* W,                                      \
         F* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ, \
         F* WORK, int64_t LWORK,                                    \
         detail::real_t<F>* RWORK, int64_t LRWORK,                  \
         scalapack_int* IWORK, int64_t LIWORK ) {                   \
                                                                    \
  auto _N      = detail::to_scalapack_int( N  );                    \
  auto _IA     = detail::to_scalapack_int( IA );                    \
  auto _JA     = detail::to_scalapack_int( JA );                    \
  auto _IZ     = detail::to_scalapack_int( IZ );                    \
  auto _JZ     = detail::to_scalapack_int( JZ );                    \
  auto _LWORK  = detail::to_scalapack_int( LWORK );                 \
  auto _LRWORK = detail::to_scalapack_int( LRWORK );                \
  auto _LIWORK = detail::to_scalapack_int( LIWORK );                \
                                                                    \
  auto _DESCA = detail::to_scalapack_int( DESCA );                  \
  auto _DESCZ = detail::to_scalapack_int( DESCZ );                  \
                                                                    \
  scalapack_int INFO;                                               \
  FNAME( JOBZ, UPLO, &_N, A, &_IA, &_JA, _DESCA.data(),             \
         W, Z, &_IZ, &_JZ, _DESCZ.data(), WORK, &_LWORK,            \
         RWORK, &_LRWORK, IWORK, &_LIWORK, &INFO );                 \
  return INFO;                                                      \
                                                                    \
}

pheevd_impl( scomplex, pcheevd_ );
pheevd_impl( dcomplex, pzheevd_ );


}

