/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/pblas/axpy.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psaxpy_( const scalapack_int* N, const float* ALPHA, const float* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX,
              float* Y, const scalapack_int* IY, const scalapack_int* JY,
              const scalapack_int* DESCY, const scalapack_int* INCY );
void pdaxpy_( const scalapack_int* N, const double* ALPHA, const double* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX,
              double* Y, const scalapack_int* IY, const scalapack_int* JY,
              const scalapack_int* DESCY, const scalapack_int* INCY );
void pcaxpy_( const scalapack_int* N, const scomplex* ALPHA, const scomplex* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX,
              scomplex* Y, const scalapack_int* IY, const scalapack_int* JY,
              const scalapack_int* DESCY, const scalapack_int* INCY );
void pzaxpy_( const scalapack_int* N, const dcomplex* ALPHA, const dcomplex* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX,
              dcomplex* Y, const scalapack_int* IY, const scalapack_int* JY,
              const scalapack_int* DESCY, const scalapack_int* INCY );

}


namespace scalapackpp {
namespace wrappers    {

#define paxpy_impl(F,FNAME)\
template <>                                                              \
void paxpy( int64_t N, F ALPHA,                                          \
  const F * X, int64_t IX, int64_t JX, const scalapack_desc& DESCX,      \
  int64_t INCX, F * Y, int64_t IY, int64_t JY,                           \
  const scalapack_desc& DESCY, int64_t INCY) {                           \
                                                                         \
  auto _N  = detail::to_scalapack_int( N  );                             \
  auto _IX = detail::to_scalapack_int( IX );                             \
  auto _JX = detail::to_scalapack_int( JX );                             \
  auto _INCX = detail::to_scalapack_int( INCX );                         \
  auto _IY = detail::to_scalapack_int( IY );                             \
  auto _JY = detail::to_scalapack_int( JY );                             \
  auto _INCY = detail::to_scalapack_int( INCY );                         \
                                                                         \
  auto _DESCX = detail::to_scalapack_int( DESCX );                       \
  auto _DESCY = detail::to_scalapack_int( DESCY );                       \
                                                                         \
  FNAME( &_N, &ALPHA, X, &_IX, &_JX, _DESCX.data(), &_INCX,              \
         Y, &_IY, &_JY, _DESCY.data(), &_INCY );                         \
                                                                         \
}


paxpy_impl( float,    psaxpy_ );
paxpy_impl( double,   pdaxpy_ );
paxpy_impl( scomplex, pcaxpy_ );
paxpy_impl( dcomplex, pzaxpy_ );
}
}
