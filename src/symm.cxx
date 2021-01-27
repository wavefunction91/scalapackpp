/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/pblas/symm.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pssymm_( const char* SIDE, const char* UPLO,
         const scalapack_int* M, const scalapack_int* N, 
         const float* ALPHA, 
         const float* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const float* B, const scalapack_int* IB, const scalapack_int* JB, 
         const scalapack_int* DESCB,
         const float* BETA,
         float* C, const scalapack_int* IC, const scalapack_int* JC, 
         const scalapack_int* DESCC );
void pdsymm_( const char* SIDE, const char* UPLO,
         const scalapack_int* M, const scalapack_int* N,  
         const double* ALPHA, 
         const double* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const double* B, const scalapack_int* IB, const scalapack_int* JB, 
         const scalapack_int* DESCB,
         const double* BETA,
         double* C, const scalapack_int* IC, const scalapack_int* JC, 
         const scalapack_int* DESCC );

}



namespace scalapackpp {
namespace wrappers    {

#define psymm_impl(F,FNAME)                                               \
template <>                                                               \
void                                                                      \
  psymm( const char* SIDE, const char* UPLO,                              \
         int64_t M, int64_t N,  F ALPHA,                                  \
         const F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, \
         const F* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB, \
         F BETA,                                                          \
         F* C, int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {    \
                                                                          \
  auto _M     = detail::to_scalapack_int( M  );                           \
  auto _N     = detail::to_scalapack_int( N  );                           \
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
  FNAME( SIDE, UPLO, &_M, &_N, &ALPHA,                                    \
         A, &_IA, &_JA, _DESCA.data(),                                    \
         B, &_IB, &_JB, _DESCB.data(), &BETA,                             \
         C, &_IC, &_JC, _DESCC.data() );                                  \
                                                                          \
}

psymm_impl( float,    pssymm_ );
psymm_impl( double,   pdsymm_ );

}
}

