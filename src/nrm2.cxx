/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/pblas/nrm2.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psnrm2_( const scalapack_int* N, float* NRM, const float* X, const scalapack_int* IX,
              const scalapack_int* JX, const scalapack_int* DESCX, 
              const scalapack_int* INCX);
void pdnrm2_( const scalapack_int* N,  double* NRM, const double* X, const scalapack_int* IX,
              const scalapack_int* JX, const scalapack_int* DESCX, 
              const scalapack_int* INCX);
void pscnrm2_( const scalapack_int* N, float* NRM, const scomplex* X, const scalapack_int* IX,
               const scalapack_int* JX, const scalapack_int* DESCX, 
               const scalapack_int* INCX);
void pdznrm2_( const scalapack_int* N,  double* NRM, const dcomplex* X, const scalapack_int* IX,
               const scalapack_int* JX, const scalapack_int* DESCX, 
               const scalapack_int* INCX);
}


namespace scalapackpp {
namespace wrappers    {

#define pnrm2_impl(BNAME, DTYPE, F, FNAME)\
template <>                                                                             \
void BNAME( int64_t N, DTYPE* NRM,                                                      \
  const F * X, int64_t IX, int64_t JX, const scalapack_desc& DESCX, int64_t INCX) {     \
                                                                                        \
  auto _N  = detail::to_scalapack_int( N  );                                            \
  auto _IX = detail::to_scalapack_int( IX );                                            \
  auto _JX = detail::to_scalapack_int( JX );                                            \
  auto _INCX = detail::to_scalapack_int( INCX );                                        \
                                                                                        \
  auto _DESCX = detail::to_scalapack_int( DESCX );                                      \
                                                                                        \
  FNAME( &_N, NRM, X, &_IX, &_JX, _DESCX.data(), &_INCX );                              \
                                                                                        \
}


pnrm2_impl( pnrm2, float,   float,    psnrm2_  );
pnrm2_impl( pnrm2, double,  double,   pdnrm2_  );
pnrm2_impl( pnrm2, float,   scomplex, pscnrm2_ );
pnrm2_impl( pnrm2, double,  dcomplex, pdznrm2_ );

}
}
