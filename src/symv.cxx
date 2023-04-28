/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/pblas/symv.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pssymv_( const char* UPLO,
         const scalapack_int* N, 
         const float* ALPHA, 
         const float* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const float* X, const scalapack_int* IX, const scalapack_int* JX, 
         const scalapack_int* DESCX, const scalapack_int* INCX,
         const float* BETA,
         float* Y, const scalapack_int* IY, const scalapack_int* JY, 
         const scalapack_int* DESCY, const scalapack_int* INCY );
void pdsymv_( const char* UPLO,
         const scalapack_int* N, 
         const double* ALPHA, 
         const double* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const double* X, const scalapack_int* IX, const scalapack_int* JX, 
         const scalapack_int* DESCX, const scalapack_int* INCX,
         const double* BETA,
         double* Y, const scalapack_int* IY, const scalapack_int* JY, 
         const scalapack_int* DESCY, const scalapack_int* INCY );
#if 0
void pcsymv_( const char* UPLO,
         const scalapack_int* N, 
         const scomplex* ALPHA, 
         const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const scomplex* X, const scalapack_int* IX, const scalapack_int* JX, 
         const scalapack_int* DESCX, const scalapack_int* INCX,
         const scomplex* BETA,
         scomplex* Y, const scalapack_int* IY, const scalapack_int* JY, 
         const scalapack_int* DESCY, const scalapack_int* INCY );
void pzsymv_( const char* UPLO,
         const scalapack_int* N, 
         const dcomplex* ALPHA, 
         const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const dcomplex* X, const scalapack_int* IX, const scalapack_int* JX, 
         const scalapack_int* DESCX, const scalapack_int* INCX,
         const dcomplex* BETA,
         dcomplex* Y, const scalapack_int* IY, const scalapack_int* JY, 
         const scalapack_int* DESCY, const scalapack_int* INCY );
#endif

}



namespace scalapackpp {
namespace wrappers    {

#define psymv_impl(F,FNAME)                                               \
template <>                                                               \
void                                                                      \
  psymv( const char* UPLO,                                                \
         int64_t N, F ALPHA,                                              \
         const F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, \
         const F* X, int64_t IX, int64_t JX, const scalapack_desc& DESCX, \
         int64_t INCX, F BETA,                                            \
         F* Y, int64_t IY, int64_t JY, const scalapack_desc& DESCY,       \
         int64_t INCY ) {    \
                                                                          \
  auto _N     = detail::to_scalapack_int( N  );                           \
  auto _IA    = detail::to_scalapack_int( IA );                           \
  auto _JA    = detail::to_scalapack_int( JA );                           \
  auto _IX    = detail::to_scalapack_int( IX );                           \
  auto _JX    = detail::to_scalapack_int( JX );                           \
  auto _IY    = detail::to_scalapack_int( IY );                           \
  auto _JY    = detail::to_scalapack_int( JY );                           \
                                                                          \
  auto _INCX    = detail::to_scalapack_int( INCX );                       \
  auto _INCY    = detail::to_scalapack_int( INCY );                       \
                                                                          \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
  auto _DESCX = detail::to_scalapack_int( DESCX );                        \
  auto _DESCY = detail::to_scalapack_int( DESCY );                        \
                                                                          \
  FNAME( UPLO, &_N, &ALPHA,                                               \
         A, &_IA, &_JA, _DESCA.data(),                                    \
         X, &_IX, &_JX, _DESCX.data(), &_INCX, &BETA,                     \
         Y, &_IY, &_JY, _DESCY.data(), &_INCY );                          \
                                                                          \
}

psymv_impl( float,    pssymv_ );
psymv_impl( double,   pdsymv_ );
//psymv_impl( scomplex, pcsymv_ );
//psymv_impl( dcomplex, pzsymv_ );

}
}

