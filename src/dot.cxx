/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/pblas/dot.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psdot_( const scalapack_int* N, float* DOT, const float* X, const scalapack_int* IX,
             const scalapack_int* JX, const scalapack_int* DESCX, 
             const scalapack_int* INCX, const float* Y, const scalapack_int* IY,
             const scalapack_int* JY, const scalapack_int* DESCY, 
             const scalapack_int* INCY );
void pddot_( const scalapack_int* N,  double* DOT, const double* X, const scalapack_int* IX,
             const scalapack_int* JX, const scalapack_int* DESCX, 
             const scalapack_int* INCX, const double* Y, const scalapack_int* IY,
             const scalapack_int* JY, const scalapack_int* DESCY, 
             const scalapack_int* INCY );

void pcdotu_( const scalapack_int* N, scomplex* DOT, const scomplex* X, const scalapack_int* IX,
              const scalapack_int* JX, const scalapack_int* DESCX, 
              const scalapack_int* INCX, const scomplex* Y, const scalapack_int* IY,
              const scalapack_int* JY, const scalapack_int* DESCY, 
              const scalapack_int* INCY );
void pzdotu_( const scalapack_int* N, dcomplex* DOT, const dcomplex* X, const scalapack_int* IX,
              const scalapack_int* JX, const scalapack_int* DESCX, 
              const scalapack_int* INCX, const dcomplex* Y, const scalapack_int* IY,
              const scalapack_int* JY, const scalapack_int* DESCY, 
              const scalapack_int* INCY );

void pcdotc_( const scalapack_int* N, scomplex* DOT, const scomplex* X, const scalapack_int* IX,
              const scalapack_int* JX, const scalapack_int* DESCX, 
              const scalapack_int* INCX, const scomplex* Y, const scalapack_int* IY,
              const scalapack_int* JY, const scalapack_int* DESCY, 
              const scalapack_int* INCY );
void pzdotc_( const scalapack_int* N, dcomplex* DOT, const dcomplex* X, const scalapack_int* IX,
              const scalapack_int* JX, const scalapack_int* DESCX, 
              const scalapack_int* INCX, const dcomplex* Y, const scalapack_int* IY,
              const scalapack_int* JY, const scalapack_int* DESCY, 
              const scalapack_int* INCY );

}


namespace scalapackpp {
namespace wrappers    {

#define pdot_impl(BNAME, DTYPE, F, FNAME)\
template <>                                                                             \
void BNAME( int64_t N, DTYPE* DOT,                                                      \
  const F * X, int64_t IX, int64_t JX, const scalapack_desc& DESCX, int64_t INCX,       \
  const F * Y, int64_t IY, int64_t JY, const scalapack_desc& DESCY, int64_t INCY  ) {   \
                                                                                        \
  auto _N  = detail::to_scalapack_int( N  );                                            \
  auto _IX = detail::to_scalapack_int( IX );                                            \
  auto _JX = detail::to_scalapack_int( JX );                                            \
  auto _INCX = detail::to_scalapack_int( INCX );                                        \
  auto _IY = detail::to_scalapack_int( IY );                                            \
  auto _JY = detail::to_scalapack_int( JY );                                            \
  auto _INCY = detail::to_scalapack_int( INCY );                                        \
                                                                                        \
  auto _DESCX = detail::to_scalapack_int( DESCX );                                      \
  auto _DESCY = detail::to_scalapack_int( DESCY );                                      \
                                                                                        \
  FNAME( &_N, DOT, X, &_IX, &_JX, _DESCX.data(), &_INCX, Y, &_IY, &_JY, _DESCY.data(),  \
         &_INCY );                                                                      \
                                                                                        \
}


pdot_impl( pdot,  float,     float,    psdot_ );
pdot_impl( pdot,  double,    double,   pddot_ );
pdot_impl( pdotu, scomplex,  scomplex, pcdotu_ );
pdot_impl( pdotu, dcomplex,  dcomplex, pzdotu_ );
pdot_impl( pdotc, scomplex,  scomplex, pcdotc_ );
pdot_impl( pdotc, dcomplex,  dcomplex, pzdotc_ );
}
}
