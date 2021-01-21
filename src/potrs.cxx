/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/linear_systems/potrs.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pspotrs_( const char* UPLO, const scalapack_int* N, const scalapack_int* NRHS, 
    const float* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    float* B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, scalapack_int* INFO );
void pdpotrs_( const char* UPLO, const scalapack_int* N, const scalapack_int* NRHS, 
    const double* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    double* B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, scalapack_int* INFO );
void pcpotrs_( const char* UPLO, const scalapack_int* N, const scalapack_int* NRHS, 
    const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    scomplex* B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, scalapack_int* INFO );
void pzpotrs_( const char* UPLO, const scalapack_int* N, const scalapack_int* NRHS, 
    const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    dcomplex* B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, scalapack_int* INFO );

}




namespace scalapackpp::wrappers {

#define ppotrs_impl(F,FNAME)                                               \
template <>                                                               \
int64_t                                                                   \
  ppotrs( const char* UPLO, int64_t N, int64_t NRHS,                       \
    const F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,      \
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
  FNAME( UPLO, &_N, &_NRHS, A, &_IA, &_JA, _DESCA.data(),                 \
    B, &_IB, &_JB, _DESCB.data(), &INFO );                                \
  return INFO;                                                            \
                                                                          \
}

ppotrs_impl( float,    pspotrs_ );
ppotrs_impl( double,   pdpotrs_ );
ppotrs_impl( scomplex, pcpotrs_ );
ppotrs_impl( dcomplex, pzpotrs_ );

}
