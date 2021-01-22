/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/pblas/gemm.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psgemm_( const char* TRANSA, const char* TRANSB,
         const scalapack_int* M, const scalapack_int* N, const scalapack_int* K, 
         const float* ALPHA, 
         const float* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const float* B, const scalapack_int* IB, const scalapack_int* JB, 
         const scalapack_int* DESCB,
         const float* BETA,
         float* C, const scalapack_int* IC, const scalapack_int* JC, 
         const scalapack_int* DESCC );
void pdgemm_( const char* TRANSA, const char* TRANSB,
         const scalapack_int* M, const scalapack_int* N, const scalapack_int* K, 
         const double* ALPHA, 
         const double* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const double* B, const scalapack_int* IB, const scalapack_int* JB, 
         const scalapack_int* DESCB,
         const double* BETA,
         double* C, const scalapack_int* IC, const scalapack_int* JC, 
         const scalapack_int* DESCC );
void pcgemm_( const char* TRANSA, const char* TRANSB,
         const scalapack_int* M, const scalapack_int* N, const scalapack_int* K, 
         const scomplex* ALPHA, 
         const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const scomplex* B, const scalapack_int* IB, const scalapack_int* JB, 
         const scalapack_int* DESCB,
         const scomplex* BETA,
         scomplex* C, const scalapack_int* IC, const scalapack_int* JC, 
         const scalapack_int* DESCC );
void pzgemm_( const char* TRANSA, const char* TRANSB,
         const scalapack_int* M, const scalapack_int* N, const scalapack_int* K, 
         const dcomplex* ALPHA, 
         const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const dcomplex* B, const scalapack_int* IB, const scalapack_int* JB, 
         const scalapack_int* DESCB,
         const dcomplex* BETA,
         dcomplex* C, const scalapack_int* IC, const scalapack_int* JC, 
         const scalapack_int* DESCC );

}



namespace scalapackpp {
namespace wrappers    {

#define pgemm_impl(F,FNAME)                                               \
template <>                                                               \
void                                                                      \
  pgemm( const char* TRANSA, const char* TRANSB,                          \
         int64_t M, int64_t N, int64_t K, F ALPHA,                        \
         const F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, \
         const F* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB, \
         F BETA,                                                          \
         F* C, int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {    \
                                                                          \
  auto _M     = detail::to_scalapack_int( M  );                           \
  auto _N     = detail::to_scalapack_int( N  );                           \
  auto _K     = detail::to_scalapack_int( K  );                           \
  auto _IA    = detail::to_scalapack_int( IA );                           \
  auto _JA    = detail::to_scalapack_int( JA );                           \
  auto _IB    = detail::to_scalapack_int( IB );                           \
  auto _JB    = detail::to_scalapack_int( JB );                           \
  auto _IC    = detail::to_scalapack_int( IC );                           \
  auto _JC    = detail::to_scalapack_int( JC );                           \
                                                                          \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
  auto _DESCB = detail::to_scalapack_int( DESCB );                        \
  auto _DESCC = detail::to_scalapack_int( DESCC );                        \
                                                                          \
  FNAME( TRANSA, TRANSB, &_M, &_N, &_K, &ALPHA,                           \
         A, &_IA, &_JA, _DESCA.data(),                                    \
         B, &_IB, &_JB, _DESCB.data(), &BETA,                             \
         C, &_IC, &_JC, _DESCC.data() );                                  \
                                                                          \
}

pgemm_impl( float,    psgemm_ );
pgemm_impl( double,   pdgemm_ );
pgemm_impl( scomplex, pcgemm_ );
pgemm_impl( dcomplex, pzgemm_ );

}
}

