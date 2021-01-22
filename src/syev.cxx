/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/eigenvalue_problem/syev.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;

// Prototypes
extern "C" {

void pssyev_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         float* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, float* W,
         float* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, float* WORK, const scalapack_int *LWORK,
         scalapack_int* INFO );

void pdsyev_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         double* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, double* W,
         double* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, double* WORK, const scalapack_int *LWORK,
         scalapack_int* INFO );

void pssyevd_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         float* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, float* W,
         float* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, float* WORK, const scalapack_int *LWORK,
         scalapack_int* IWORK, const scalapack_int* LIWORK, scalapack_int* INFO );

void pdsyevd_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         double* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, double* W,
         double* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, double* WORK, const scalapack_int *LWORK,
         scalapack_int* IWORK, const scalapack_int* LIWORK, scalapack_int* INFO );

}


namespace scalapackpp {
namespace wrappers    {

#define psyev_impl(F, FNAME)                                        \
template <>                                                         \
int64_t                                                             \
  psyev( const char* JOBZ, const char* UPLO, int64_t N,             \
         F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, \
         F* W,                                                      \
         F* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ, \
         F* WORK, int64_t LWORK ) {                                 \
                                                                    \
  auto _N     = detail::to_scalapack_int( N  );                     \
  auto _IA    = detail::to_scalapack_int( IA );                     \
  auto _JA    = detail::to_scalapack_int( JA );                     \
  auto _IZ    = detail::to_scalapack_int( IZ );                     \
  auto _JZ    = detail::to_scalapack_int( JZ );                     \
  auto _LWORK = detail::to_scalapack_int( LWORK );                  \
                                                                    \
  auto _DESCA = detail::to_scalapack_int( DESCA );                  \
  auto _DESCZ = detail::to_scalapack_int( DESCZ );                  \
                                                                    \
  scalapack_int INFO;                                               \
  FNAME( JOBZ, UPLO, &_N, A, &_IA, &_JA, _DESCA.data(),             \
         W, Z, &_IZ, &_JZ, _DESCZ.data(), WORK, &_LWORK, &INFO );   \
  return INFO;                                                      \
                                                                    \
}

psyev_impl( float,  pssyev_ );
psyev_impl( double, pdsyev_ );

#define psyevd_impl(F, FNAME)                                       \
template <>                                                         \
int64_t                                                             \
  psyevd( const char* JOBZ, const char* UPLO, int64_t N,            \
         F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, \
         F* W,                                                      \
         F* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ, \
         F* WORK, int64_t LWORK,                                    \
         scalapack_int* IWORK, int64_t LIWORK ) {                   \
                                                                    \
  auto _N      = detail::to_scalapack_int( N  );                    \
  auto _IA     = detail::to_scalapack_int( IA );                    \
  auto _JA     = detail::to_scalapack_int( JA );                    \
  auto _IZ     = detail::to_scalapack_int( IZ );                    \
  auto _JZ     = detail::to_scalapack_int( JZ );                    \
  auto _LWORK  = detail::to_scalapack_int( LWORK );                 \
  auto _LIWORK = detail::to_scalapack_int( LIWORK );                \
                                                                    \
  auto _DESCA = detail::to_scalapack_int( DESCA );                  \
  auto _DESCZ = detail::to_scalapack_int( DESCZ );                  \
                                                                    \
  scalapack_int INFO;                                               \
  FNAME( JOBZ, UPLO, &_N, A, &_IA, &_JA, _DESCA.data(),             \
         W, Z, &_IZ, &_JZ, _DESCZ.data(), WORK, &_LWORK,            \
         IWORK, &_LIWORK, &INFO );                                  \
  return INFO;                                                      \
                                                                    \
}

psyevd_impl( float,  pssyevd_ );
psyevd_impl( double, pdsyevd_ );

}
}

