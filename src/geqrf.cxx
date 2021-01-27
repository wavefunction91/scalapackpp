/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/factorizations/geqrf.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psgeqrf_( const scalapack_int* M, const scalapack_int* N, float* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  float* TAU, float* WORK, const scalapack_int* LWORK, scalapack_int* INFO );
void pdgeqrf_( const scalapack_int* M, const scalapack_int* N, double* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  double* TAU, double* WORK, const scalapack_int* LWORK, scalapack_int* INFO );
void pcgeqrf_( const scalapack_int* M, const scalapack_int* N, scomplex* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scomplex* TAU, scomplex* WORK, const scalapack_int* LWORK, scalapack_int* INFO );
void pzgeqrf_( const scalapack_int* M, const scalapack_int* N, dcomplex* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  dcomplex* TAU, dcomplex* WORK, const scalapack_int* LWORK, scalapack_int* INFO );

}


namespace scalapackpp {
namespace wrappers    {

#define pgeqrf_impl(F,FNAME)                                              \
template <>                                                               \
int64_t                                                                   \
  pgeqrf( int64_t M, int64_t N, F* A, int64_t IA, int64_t JA,             \
          const scalapack_desc& DESCA, F* TAU, F* WORK,                   \
          int64_t LWORK ) {                                               \
                                                                          \
  auto _M     = detail::to_scalapack_int( M    );                         \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _LWORK = detail::to_scalapack_int( LWORK );                        \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( &_M, &_N, A, &_IA, &_JA, _DESCA.data(), TAU, WORK, &_LWORK,      \
         &INFO );                                                         \
  return INFO;                                                            \
                                                                          \
}

pgeqrf_impl( float,    psgeqrf_ );
pgeqrf_impl( double,   pdgeqrf_ );
pgeqrf_impl( scomplex, pcgeqrf_ );
pgeqrf_impl( dcomplex, pzgeqrf_ );


}
}

