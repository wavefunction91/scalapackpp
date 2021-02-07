/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/matrix_inverse/potri.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pspotri_( const char* UPLO, const scalapack_int* N, float* A, const scalapack_int* IA, 
          const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* INFO );
void pdpotri_( const char* UPLO, const scalapack_int* N, double* A, const scalapack_int* IA, 
          const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* INFO );
void pcpotri_( const char* UPLO, const scalapack_int* N, scomplex* A, const scalapack_int* IA, 
          const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* INFO );
void pzpotri_( const char* UPLO, const scalapack_int* N, dcomplex* A, const scalapack_int* IA, 
          const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* INFO );

}

namespace scalapackpp {
namespace wrappers    {

#define ppotri_impl(F,FNAME)                                              \
template <>                                                               \
int64_t                                                                   \
  ppotri( const char* UPLO, int64_t N, F* A, int64_t IA, int64_t JA,      \
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

ppotri_impl( float,    pspotri_ );
ppotri_impl( double,   pdpotri_ );
ppotri_impl( scomplex, pcpotri_ );
ppotri_impl( dcomplex, pzpotri_ );

}
}
