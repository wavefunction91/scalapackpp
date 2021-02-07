/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/gemr2d.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;


extern "C" {


void psgemr2d_( const scalapack_int* M, const scalapack_int* N, 
    const float* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    float * B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, const scalapack_int* ICONTEXT );
void pdgemr2d_( const scalapack_int* M, const scalapack_int* N, 
    const double* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    double * B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, const scalapack_int* ICONTEXT );
void pcgemr2d_( const scalapack_int* M, const scalapack_int* N, 
    const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    scomplex * B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, const scalapack_int* ICONTEXT );
void pzgemr2d_( const scalapack_int* M, const scalapack_int* N, 
    const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    dcomplex * B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, const scalapack_int* ICONTEXT );

}


namespace scalapackpp {
namespace wrappers    {

#define pgemr2d_impl(F, FNAME)                                            \
template <>                                                               \
void pgemr2d( int64_t M, int64_t N,                                       \
    const F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,      \
    F* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB,            \
    int64_t ICONTEXT ) {                                                  \
                                                                          \
  auto _M     = detail::to_scalapack_int( M  );                           \
  auto _N     = detail::to_scalapack_int( N  );                           \
  auto _IA    = detail::to_scalapack_int( IA );                           \
  auto _JA    = detail::to_scalapack_int( JA );                           \
  auto _IB    = detail::to_scalapack_int( IB );                           \
  auto _JB    = detail::to_scalapack_int( JB );                           \
                                                                          \
  auto _ICONTEXT = detail::to_scalapack_int( ICONTEXT );                  \
                                                                          \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
  auto _DESCB = detail::to_scalapack_int( DESCB );                        \
                                                                          \
  FNAME( &_M, &_N, A, &_IA, &_JA, _DESCA.data(),                          \
                   B, &_IB, &_JB, _DESCB.data(), &_ICONTEXT );            \
                                                                          \
}

pgemr2d_impl( float,    psgemr2d_ );
pgemr2d_impl( double,   pdgemr2d_ );
pgemr2d_impl( scomplex, pcgemr2d_ );
pgemr2d_impl( dcomplex, pzgemr2d_ );

}
}
