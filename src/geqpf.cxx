/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/factorizations/geqpf.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psgeqpf_( const scalapack_int* M, const scalapack_int* N, float* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, float* TAU, float* WORK, const scalapack_int* LWORK, 
  scalapack_int* INFO );
void pdgeqpf_( const scalapack_int* M, const scalapack_int* N, double* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, double* TAU, double* WORK, const scalapack_int* LWORK, 
  scalapack_int* INFO );
void pcgeqpf_( const scalapack_int* M, const scalapack_int* N, scomplex* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, scomplex* TAU, scomplex* WORK, const scalapack_int* LWORK, 
  float* RWORK, const scalapack_int* LRWORK, scalapack_int* INFO );
void pzgeqpf_( const scalapack_int* M, const scalapack_int* N, dcomplex* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, dcomplex* TAU, dcomplex* WORK, const scalapack_int* LWORK, 
  double* RWORK, const scalapack_int* LRWORK, scalapack_int* INFO );

}


namespace scalapackpp {
namespace wrappers    {

#define pgeqpf_real_impl(F,FNAME)                                         \
template <>                                                               \
int64_t                                                                   \
  pgeqpf( int64_t M, int64_t N, F* A, int64_t IA, int64_t JA,             \
          const scalapack_desc& DESCA, scalapack_int* IPIV, F* TAU,       \
          F* WORK, int64_t LWORK ) {                                      \
                                                                          \
  auto _M     = detail::to_scalapack_int( M    );                         \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _LWORK = detail::to_scalapack_int( LWORK );                        \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( &_M, &_N, A, &_IA, &_JA, _DESCA.data(), IPIV, TAU,               \
         WORK, &_LWORK, &INFO );                                          \
  return INFO;                                                            \
                                                                          \
}

pgeqpf_real_impl( float,    psgeqpf_ );
pgeqpf_real_impl( double,   pdgeqpf_ );

#define pgeqpf_complex_impl(F,FNAME)                                      \
template <>                                                               \
int64_t                                                                   \
  pgeqpf( int64_t M, int64_t N, F* A, int64_t IA, int64_t JA,             \
          const scalapack_desc& DESCA, scalapack_int* IPIV, F* TAU,       \
          F* WORK, int64_t LWORK, detail::real_t<F>* RWORK,               \
          int64_t LRWORK) {                                               \
                                                                          \
  auto _M     = detail::to_scalapack_int( M    );                         \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _LWORK = detail::to_scalapack_int( LWORK );                        \
  auto _LRWORK = detail::to_scalapack_int( LRWORK );                      \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( &_M, &_N, A, &_IA, &_JA, _DESCA.data(), IPIV, TAU,               \
         WORK, &_LWORK, RWORK, &_LRWORK, &INFO );                         \
  return INFO;                                                            \
                                                                          \
}

pgeqpf_complex_impl( scomplex, pcgeqpf_ );
pgeqpf_complex_impl( dcomplex, pzgeqpf_ );


}
}

