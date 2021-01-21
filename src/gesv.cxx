/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/linear_systems/gesv.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psgesv_( const scalapack_int* N, const scalapack_int* NRHS, 
    float* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    float* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pdgesv_( const scalapack_int* N, const scalapack_int* NRHS, 
    double* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    double* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pcgesv_( const scalapack_int* N, const scalapack_int* NRHS, 
    scomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    scomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pzgesv_( const scalapack_int* N, const scalapack_int* NRHS, 
    dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    dcomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );

}


namespace scalapackpp::wrappers {

#define pgesv_impl(F,FNAME)                                               \
template <>                                                               \
int64_t                                                                   \
  pgesv( int64_t N, int64_t NRHS,                                         \
    F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,            \
    const scalapack_int* IPIV,                                            \
    F* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {         \
                                                                          \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _NRHS  = detail::to_scalapack_int( NRHS );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _IB    = detail::to_scalapack_int( IB   );                         \
  auto _JB    = detail::to_scalapack_int( JB   );                         \
                                                                          \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
  auto _DESCB = detail::to_scalapack_int( DESCB );                        \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( &_N, &_NRHS, A, &_IA, &_JA, _DESCA.data(), IPIV,                 \
    B, &_IB, &_JB, _DESCB.data(), &INFO );                                \
  return INFO;                                                            \
                                                                          \
}

pgesv_impl( float,    psgesv_ );
pgesv_impl( double,   pdgesv_ );
pgesv_impl( scomplex, pcgesv_ );
pgesv_impl( dcomplex, pzgesv_ );
          
}
