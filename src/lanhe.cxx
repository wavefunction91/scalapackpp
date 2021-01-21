/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/matrix_norm/lanhe.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

float pclanhe_( const char* NORM, const char* UPLO, const scalapack_int* N,
                const scomplex* A, const scalapack_int* IA, const scalapack_int* JA,
                const scalapack_int* DESCA, float* WORK );
double pzlanhe_( const char* NORM, const char* UPLO, const scalapack_int* N,
                 const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA,
                 const scalapack_int* DESCA, double* WORK );

}

namespace scalapackpp::wrappers {

#define planhe_impl(F,FNAME)                                              \
template <>                                                               \
detail::real_t<F> planhe( const char* NORM, const char* UPLO, int64_t N,  \
          const F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,\
          detail::real_t<F>* WORK ) {                                     \
                                                                          \
  auto _N     = detail::to_scalapack_int( N  );                           \
  auto _IA    = detail::to_scalapack_int( IA );                           \
  auto _JA    = detail::to_scalapack_int( JA );                           \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
                                                                          \
  return FNAME( NORM, UPLO, &_N, A, &_IA, &_JA, _DESCA.data(), WORK );    \
                                                                          \
}

planhe_impl( scomplex, pclanhe_ );
planhe_impl( dcomplex, pzlanhe_ );

}
