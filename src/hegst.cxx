/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/eigenvalue_problem/hegst.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pchegst_( const scalapack_int* IBTYPE, const char* UPLO, const scalapack_int* N,
               scomplex* A, const scalapack_int* IA, const scalapack_int* JA,
               const scalapack_int* DESCA,
               const scomplex* B, const scalapack_int* IB, const scalapack_int* JB,
               const scalapack_int* DESCB,
               float* SCALE, scalapack_int* INFO );
void pzhegst_( const scalapack_int* IBTYPE, const char* UPLO, const scalapack_int* N,
               dcomplex* A, const scalapack_int* IA, const scalapack_int* JA,
               const scalapack_int* DESCA,
               const dcomplex* B, const scalapack_int* IB, const scalapack_int* JB,
               const scalapack_int* DESCB,
               double* SCALE, scalapack_int* INFO );

}

namespace scalapackpp {
namespace wrappers    {

#define phegst_impl(F, FNAME)                                              \
template <>                                                                \
int64_t                                                                    \
  phegst( int64_t IBTYPE, const char* UPLO, int64_t N,                     \
                F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, \
          const F* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB, \
                detail::real_t<F>* SCALE ) {                               \
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

phegst_impl( scomplex, pchegst_ );
phegst_impl( dcomplex, pzhegst_ );

}
}
