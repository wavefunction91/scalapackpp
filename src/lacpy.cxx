/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/lacpy.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pslacpy_( const char*, const scalapack_int* M, const scalapack_int* N,
               const float* A, const scalapack_int* IA, const scalapack_int* JA,
               const scalapack_int* DESCA, float* B, const scalapack_int* IB,
               const scalapack_int* JB, const scalapack_int* DESCB );

void pdlacpy_( const char*, const scalapack_int* M, const scalapack_int* N,
               const double* A, const scalapack_int* IA, const scalapack_int* JA,
               const scalapack_int* DESCA, double* B, const scalapack_int* IB,
               const scalapack_int* JB, const scalapack_int* DESCB );

void pclacpy_( const char*, const scalapack_int* M, const scalapack_int* N,
               const scomplex* A, const scalapack_int* IA, const scalapack_int* JA,
               const scalapack_int* DESCA, scomplex* B, const scalapack_int* IB,
               const scalapack_int* JB, const scalapack_int* DESCB );

void pzlacpy_( const char*, const scalapack_int* M, const scalapack_int* N,
               const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA,
               const scalapack_int* DESCA, dcomplex* B, const scalapack_int* IB,
               const scalapack_int* JB, const scalapack_int* DESCB );

}

namespace scalapackpp {
namespace wrappers    {

#define placpy_impl(F, FNAME)                                                 \
template <>                                                                   \
void placpy( const char* UPLO, int64_t M, int64_t N, const F* A, int64_t IA,  \
             int64_t JA, const scalapack_desc& DESCA, F* B, int64_t IB,       \
             int64_t JB, const scalapack_desc& DESCB ) {                      \
                                                                              \
  auto _M     = detail::to_scalapack_int( M  );                               \
  auto _N     = detail::to_scalapack_int( N  );                               \
  auto _IA    = detail::to_scalapack_int( IA );                               \
  auto _JA    = detail::to_scalapack_int( JA );                               \
  auto _IB    = detail::to_scalapack_int( IB );                               \
  auto _JB    = detail::to_scalapack_int( JB );                               \
                                                                              \
  auto _DESCA = detail::to_scalapack_int( DESCA );                            \
  auto _DESCB = detail::to_scalapack_int( DESCB );                            \
                                                                              \
  FNAME( UPLO, &_M, &_N, A, &_IA, &_JA, _DESCA.data(),                        \
         B, &_IB, &_JB, _DESCB.data() );                                      \
}

placpy_impl( float,    pslacpy_ );
placpy_impl( double,   pdlacpy_ );
placpy_impl( scomplex, pclacpy_ );
placpy_impl( dcomplex, pzlacpy_ );

}
}
