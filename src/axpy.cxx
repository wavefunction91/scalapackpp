/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/pblas/axpy.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psaxpy_( const scalapack_int* N, const float* ALPHA, const float* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, float* Y, const scalapack_int* IY, 
              const scalapack_int* JY, const scalapack_int* DESCY );
void pdaxpy_( const scalapack_int* N, const double* ALPHA, const double* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, double* Y, const scalapack_int* IY, 
              const scalapack_int* JY, const scalapack_int* DESCY );
void pcaxpy_( const scalapack_int* N, const scomplex* ALPHA, const scomplex* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, scomplex* Y, const scalapack_int* IY, 
              const scalapack_int* JY, const scalapack_int* DESCY );
void pzaxpy_( const scalapack_int* N, const dcomplex* ALPHA, const dcomplex* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, dcomplex* Y, const scalapack_int* IY, 
              const scalapack_int* JY, const scalapack_int* DESCY );

}


namespace scalapackpp::wrappers {

#define paxpy_impl(F,FNAME)\
template <>                                                              \
void paxpy( int64_t N, float ALPHA,                                      \
    const F * X, int64_t IX, int64_t JX, const scalapack_desc& DESCX,    \
          F * Y, int64_t IY, int64_t JY, const scalapack_desc& DESCY ) { \
                                                                         \
  auto _N  = detail::to_scalapack_int( N  );                             \
  auto _IX = detail::to_scalapack_int( IX );                             \
  auto _JX = detail::to_scalapack_int( JX );                             \
  auto _IY = detail::to_scalapack_int( IY );                             \
  auto _JY = detail::to_scalapack_int( JY );                             \
                                                                         \
  auto _DESCX = detail::to_scalapack_int( DESCX );                       \
  auto _DESCY = detail::to_scalapack_int( DESCY );                       \
                                                                         \
  FNAME( &_N, &ALPHA, X, &_IX, &_JX, _DESCX.data(),                      \
                      Y, &_IY, &_JY, _DESCY.data() );                    \
                                                                         \
}

paxpy_impl( float,    psaxpy_ );
paxpy_impl( double,   pdaxpy_ );
paxpy_impl( scomplex, pcaxpy_ );
paxpy_impl( dcomplex, pzaxpy_ );



}
