/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/factorizations/getrf.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psgetrf_( const scalapack_int* M, const scalapack_int* N, float* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, scalapack_int* INFO );
void pdgetrf_( const scalapack_int* M, const scalapack_int* N, double* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, scalapack_int* INFO );
void pcgetrf_( const scalapack_int* M, const scalapack_int* N, scomplex* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, scalapack_int* INFO );
void pzgetrf_( const scalapack_int* M, const scalapack_int* N, dcomplex* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, scalapack_int* INFO );

}


namespace scalapackpp::wrappers {

#define pgetrf_impl(F,FNAME)                                              \
template <>                                                               \
int64_t                                                                   \
  pgetrf( int64_t M, int64_t N, F* A, int64_t IA, int64_t JA,             \
          const scalapack_desc& DESCA, scalapack_int* IPIV ) {            \
                                                                          \
  auto _M     = detail::to_scalapack_int( M    );                         \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( &_M, &_N, A, &_IA, &_JA, _DESCA.data(), IPIV, &INFO );           \
  return INFO;                                                            \
                                                                          \
}

pgetrf_impl( float,    psgetrf_ );
pgetrf_impl( double,   pdgetrf_ );
pgetrf_impl( scomplex, pcgetrf_ );
pgetrf_impl( dcomplex, pzgetrf_ );


}
