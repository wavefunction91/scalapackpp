/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/eigenvalue_problem/sygst.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;

// Prototypes
extern "C" {

void pssygst_( const scalapack_int* IBTYPE, const char* UPLO, const scalapack_int* N,
               float* A, const scalapack_int* IA, const scalapack_int* JA,
               const scalapack_int* DESCA,
               const float* B, const scalapack_int* IB, const scalapack_int* JB,
               const scalapack_int* DESCB,
               float* SCALE, scalapack_int* INFO );
void pdsygst_( const scalapack_int* IBTYPE, const char* UPLO, const scalapack_int* N,
               double* A, const scalapack_int* IA, const scalapack_int* JA,
               const scalapack_int* DESCA,
               const double* B, const scalapack_int* IB, const scalapack_int* JB,
               const scalapack_int* DESCB,
               double* SCALE, scalapack_int* INFO );

}

namespace scalapackpp {
namespace wrappers    {

#define psygst_impl(F, FNAME)                                              \
template <>                                                                \
int64_t                                                                    \
  psygst( int64_t IBTYPE, const char* UPLO, int64_t N,                     \
                F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, \
          const F* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB, \
                F* SCALE ) {                                               \
                                                                           \
  auto _IBTYPE     = detail::to_scalapack_int( IBTYPE  );                  \
                                                                           \
  auto _N     = detail::to_scalapack_int( N  );                            \
  auto _IA    = detail::to_scalapack_int( IA );                            \
  auto _JA    = detail::to_scalapack_int( JA );                            \
  auto _IB    = detail::to_scalapack_int( IB );                            \
  auto _JB    = detail::to_scalapack_int( JB );                            \
                                                                           \
  auto _DESCA = detail::to_scalapack_int( DESCA );                         \
  auto _DESCB = detail::to_scalapack_int( DESCB );                         \
                                                                           \
  scalapack_int INFO;                                                      \
  FNAME( &_IBTYPE, UPLO, &_N, A, &_IA, &_JA, _DESCA.data(), B, &_IB, &_JB, \
         _DESCB.data(), SCALE, &INFO );                                    \
  return INFO;                                                             \
                                                                           \
}

psygst_impl( float,  pssygst_ );
psygst_impl( double, pdsygst_ );

}
}
