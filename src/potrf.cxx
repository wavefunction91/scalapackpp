/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/factorizations/potrf.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {


void pspotrf_( const char* UPLO, const scalapack_int* N, 
               float* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* );
void pdpotrf_( const char* UPLO, const scalapack_int* N, 
               double* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* );
void pcpotrf_( const char* UPLO, const scalapack_int* N, 
               scomplex* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* );
void pzpotrf_( const char* UPLO, const scalapack_int* N, 
               dcomplex* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* );

}

namespace scalapackpp {
namespace wrappers    {

#define ppotrf_impl(F,FNAME)                                              \
template <>                                                               \
int64_t                                                                   \
  ppotrf( const char* UPLO, int64_t N, F* A, int64_t IA, int64_t JA,      \
          const scalapack_desc& DESCA ) {                                 \
                                                                          \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( UPLO, &_N, A, &_IA, &_JA, _DESCA.data(), &INFO );                \
  return INFO;                                                            \
                                                                          \
}

ppotrf_impl( float,    pspotrf_ );
ppotrf_impl( double,   pdpotrf_ );
ppotrf_impl( scomplex, pcpotrf_ );
ppotrf_impl( dcomplex, pzpotrf_ );

}
}
