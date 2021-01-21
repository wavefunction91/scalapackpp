/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/matrix_norm/lange.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;


// Prototypes
extern "C" {

float pslange_( const char* NORM, const scalapack_int* M, const scalapack_int* N,
                const float* A, const scalapack_int* IA, const scalapack_int* JA,
                const scalapack_int* DESCA, float* WORK );
double pdlange_( const char* NORM, const scalapack_int* M, const scalapack_int* N,
                 const double* A, const scalapack_int* IA, const scalapack_int* JA,
                 const scalapack_int* DESCA, double* WORK );
float pclange_( const char* NORM, const scalapack_int* M, const scalapack_int* N,
                const scomplex* A, const scalapack_int* IA, const scalapack_int* JA,
                const scalapack_int* DESCA, float* WORK );
double pzlange_( const char* NORM, const scalapack_int* M, const scalapack_int* N,
                 const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA,
                 const scalapack_int* DESCA, double* WORK );

}

namespace scalapackpp::wrappers {

#define plange_impl(F,FNAME)                                              \
template <>                                                               \
detail::real_t<F> plange( const char* NORM, int64_t M, int64_t N,         \
          const F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,\
          detail::real_t<F>* WORK ) {                                     \
                                                                          \
  auto _M     = detail::to_scalapack_int( M  );                           \
  auto _N     = detail::to_scalapack_int( N  );                           \
  auto _IA    = detail::to_scalapack_int( IA );                           \
  auto _JA    = detail::to_scalapack_int( JA );                           \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
                                                                          \
  return FNAME( NORM, &_M, &_N, A, &_IA, &_JA, _DESCA.data(), WORK );     \
                                                                          \
}

plange_impl( float,    pslange_ );
plange_impl( double,   pdlange_ );
plange_impl( scomplex, pclange_ );
plange_impl( dcomplex, pzlange_ );

}
