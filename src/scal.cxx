/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/pblas/scal.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psscal_( const scalapack_int* N, const float* ALPHA, float* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX );
void pdscal_( const scalapack_int* N, const double* ALPHA, double* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX );
void pcscal_( const scalapack_int* N, const scomplex* ALPHA, scomplex* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX );
void pzscal_( const scalapack_int* N, const dcomplex* ALPHA, dcomplex* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX );

}


namespace scalapackpp::wrappers {

#define pscal_impl(F,FNAME)\
template <>                                                              \
void pscal( int64_t N, F ALPHA,                                          \
  F * X, int64_t IX, int64_t JX, const scalapack_desc& DESCX,            \
  int64_t INCX ) {                                                       \
                                                                         \
  auto _N  = detail::to_scalapack_int( N  );                             \
  auto _IX = detail::to_scalapack_int( IX );                             \
  auto _JX = detail::to_scalapack_int( JX );                             \
  auto _INCX = detail::to_scalapack_int( INCX );                         \
                                                                         \
  auto _DESCX = detail::to_scalapack_int( DESCX );                       \
                                                                         \
  FNAME( &_N, &ALPHA, X, &_IX, &_JX, _DESCX.data(), &_INCX );            \
                                                                         \
}


pscal_impl( float,    psscal_ );
pscal_impl( double,   pdscal_ );
pscal_impl( scomplex, pcscal_ );
pscal_impl( dcomplex, pzscal_ );
}
