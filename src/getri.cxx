/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/matrix_inverse/getri.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psgetri_( const scalapack_int* N, float * A, const scalapack_int* IA, const scalapack_int* JA, 
          const scalapack_int* DESCA, const scalapack_int* IPIV, 
          float* WORK, const scalapack_int* LWORK, scalapack_int* IWORK, const scalapack_int* LIWORK, 
          scalapack_int* INFO );
void pdgetri_( const scalapack_int* N, double * A, const scalapack_int* IA, const scalapack_int* JA, 
          const scalapack_int* DESCA, const scalapack_int* IPIV, 
          double* WORK, const scalapack_int* LWORK, scalapack_int* IWORK, const scalapack_int* LIWORK, 
          scalapack_int* INFO );
void pcgetri_( const scalapack_int* N, scomplex * A, const scalapack_int* IA, const scalapack_int* JA, 
          const scalapack_int* DESCA, const scalapack_int* IPIV, 
          scomplex* WORK, const scalapack_int* LWORK, scalapack_int* IWORK, const scalapack_int* LIWORK, 
          scalapack_int* INFO );
void pzgetri_( const scalapack_int* N, dcomplex * A, const scalapack_int* IA, const scalapack_int* JA, 
          const scalapack_int* DESCA, const scalapack_int* IPIV, 
          dcomplex* WORK, const scalapack_int* LWORK, scalapack_int* IWORK, const scalapack_int* LIWORK, 
          scalapack_int* INFO );

}


namespace scalapackpp {
namespace wrappers    {

#define pgetri_impl(F,FNAME)                                              \
template <>                                                               \
int64_t                                                                   \
  pgetri( int64_t N, F* A, int64_t IA, int64_t JA,                        \
          const scalapack_desc& DESCA, const scalapack_int* IPIV,         \
          F* WORK, int64_t LWORK,                                         \
          scalapack_int* IWORK, int64_t LIWORK ) {                        \
                                                                          \
  auto _N      = detail::to_scalapack_int( N    );                        \
  auto _IA     = detail::to_scalapack_int( IA   );                        \
  auto _JA     = detail::to_scalapack_int( JA   );                        \
  auto _DESCA  = detail::to_scalapack_int( DESCA );                       \
  auto _LWORK  = detail::to_scalapack_int( LWORK );                       \
  auto _LIWORK = detail::to_scalapack_int( LIWORK );                      \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( &_N, A, &_IA, &_JA, _DESCA.data(), IPIV, WORK, &_LWORK,          \
         IWORK, &_LIWORK, &INFO );                                        \
  return INFO;                                                            \
                                                                          \
}

pgetri_impl( float,    psgetri_ );
pgetri_impl( double,   pdgetri_ );
pgetri_impl( scomplex, pcgetri_ );
pgetri_impl( dcomplex, pzgetri_ );

}
}




