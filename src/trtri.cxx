/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/matrix_inverse/trtri.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pstrtri_( const char* UPLO, const char* DIAG, const scalapack_int* N,
          float* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
          scalapack_int* INFO );
void pdtrtri_( const char* UPLO, const char* DIAG, const scalapack_int* N,
          double* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
          scalapack_int* INFO );
void pctrtri_( const char* UPLO, const char* DIAG, const scalapack_int* N,
          scomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
          scalapack_int* INFO );
void pztrtri_( const char* UPLO, const char* DIAG, const scalapack_int* N,
          dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
          scalapack_int* INFO );

}


namespace scalapackpp::wrappers {

#define ptrtri_impl(F,FNAME)                                              \
template <>                                                               \
int64_t                                                                   \
  ptrtri( const char* UPLO, const char* DIAG, int64_t N,                  \
          F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA ) {   \
                                                                          \
  auto _N     = detail::to_scalapack_int( N  );                           \
  auto _IA    = detail::to_scalapack_int( IA );                           \
  auto _JA    = detail::to_scalapack_int( JA );                           \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( UPLO, DIAG, &_N, A, &_IA, &_JA, _DESCA.data(), &INFO );          \
  return INFO;                                                            \
                                                                          \
}

ptrtri_impl( float,    pstrtri_ );
ptrtri_impl( double,   pdtrtri_ );
ptrtri_impl( scomplex, pctrtri_ );
ptrtri_impl( dcomplex, pztrtri_ );

}
