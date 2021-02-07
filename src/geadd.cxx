/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/geadd.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psgeadd_( const char* TRANS, const scalapack_int* M, const scalapack_int* N, 
               const float* ALPHA, const float* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA,
               const float* BETA, float* C, const scalapack_int* IC, 
               const scalapack_int* JC, const scalapack_int* DESCC );
void pdgeadd_( const char* TRANS, const scalapack_int* M, const scalapack_int* N, 
               const double* ALPHA, const double* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA,
               const double* BETA, double* C, const scalapack_int* IC, 
               const scalapack_int* JC, const scalapack_int* DESCC );
void pcgeadd_( const char* TRANS, const scalapack_int* M, const scalapack_int* N, 
               const scomplex* ALPHA, const scomplex* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA,
               const scomplex* BETA, scomplex* C, const scalapack_int* IC, 
               const scalapack_int* JC, const scalapack_int* DESCC );
void pzgeadd_( const char* TRANS, const scalapack_int* M, const scalapack_int* N, 
               const dcomplex* ALPHA, const dcomplex* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA,
               const dcomplex* BETA, dcomplex* C, const scalapack_int* IC, 
               const scalapack_int* JC, const scalapack_int* DESCC );

}



namespace scalapackpp {
namespace wrappers    {

#define pgeadd_impl(F,FNAME)                                              \
template <>                                                               \
void                                                                      \
  pgeadd( const char* TRANS,                                              \
         int64_t M, int64_t N,  F ALPHA,                                  \
         const F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA, \
         F BETA,                                                          \
         F* C, int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {    \
                                                                          \
  auto _M     = detail::to_scalapack_int( M  );                           \
  auto _N     = detail::to_scalapack_int( N  );                           \
  auto _IA    = detail::to_scalapack_int( IA );                           \
  auto _JA    = detail::to_scalapack_int( JA );                           \
  auto _IC    = detail::to_scalapack_int( IC );                           \
  auto _JC    = detail::to_scalapack_int( JC );                           \
                                                                          \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
  auto _DESCC = detail::to_scalapack_int( DESCC );                        \
                                                                          \
  FNAME( TRANS, &_M, &_N, &ALPHA,                                         \
         A, &_IA, &_JA, _DESCA.data(),                                    \
         &BETA,                                                           \
         C, &_IC, &_JC, _DESCC.data() );                                  \
                                                                          \
}

pgeadd_impl( float,    psgeadd_ );
pgeadd_impl( double,   pdgeadd_ );
pgeadd_impl( scomplex, pcgeadd_ );
pgeadd_impl( dcomplex, pzgeadd_ );

}
}
