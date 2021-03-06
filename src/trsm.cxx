/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/pblas/trsm.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pstrsm_( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         const scalapack_int* M, const scalapack_int* N, const float* ALPHA, 
         const float* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
         float* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB );
void pdtrsm_( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         const scalapack_int* M, const scalapack_int* N, const double* ALPHA, 
         const double* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
         double* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB );
void pctrsm_( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         const scalapack_int* M, const scalapack_int* N, const scomplex* ALPHA, 
         const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
         scomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB );
void pztrsm_( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         const scalapack_int* M, const scalapack_int* N, const dcomplex* ALPHA, 
         const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
         dcomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB );

}

namespace scalapackpp {
namespace wrappers    {

#define ptrsm_impl(F,FNAME)                                               \
template <>                                                               \
void ptrsm(                                                               \
  const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,\
  int64_t M, int64_t N, F ALPHA,                                          \
  const F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,        \
  F* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {           \
                                                                          \
  auto _M     = detail::to_scalapack_int( M  );                           \
  auto _N     = detail::to_scalapack_int( N  );                           \
  auto _IA    = detail::to_scalapack_int( IA );                           \
  auto _JA    = detail::to_scalapack_int( JA );                           \
  auto _IB    = detail::to_scalapack_int( IB );                           \
  auto _JB    = detail::to_scalapack_int( JB );                           \
                                                                          \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
  auto _DESCB = detail::to_scalapack_int( DESCB );                        \
                                                                          \
  FNAME( SIDE, UPLO, TRANS, DIAG, &_M, &_N, &ALPHA,                       \
         A, &_IA, &_JA, _DESCA.data(),                                    \
         B, &_IB, &_JB, _DESCB.data() );                                  \
                                                                          \
}

ptrsm_impl( float,    pstrsm_ );
ptrsm_impl( double,   pdtrsm_ );
ptrsm_impl( scomplex, pctrsm_ );
ptrsm_impl( dcomplex, pztrsm_ );

}
}
