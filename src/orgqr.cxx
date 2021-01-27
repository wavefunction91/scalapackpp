/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/householder/orgqr.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psorgqr_( const scalapack_int* M, const scalapack_int* N, const scalapack_int* K, float* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  const float* TAU, float* WORK, const scalapack_int* LWORK, scalapack_int* INFO );
void pdorgqr_( const scalapack_int* M, const scalapack_int* N, const scalapack_int* K, double* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  const double* TAU, double* WORK, const scalapack_int* LWORK, scalapack_int* INFO );

}


namespace scalapackpp {
namespace wrappers    {

#define porgqr_impl(F,FNAME)                                              \
template <>                                                               \
int64_t                                                                   \
  porgqr( int64_t M, int64_t N, int64_t K, F* A, int64_t IA, int64_t JA,  \
          const scalapack_desc& DESCA, const F* TAU, F* WORK,             \
          int64_t LWORK ) {                                               \
                                                                          \
  auto _M     = detail::to_scalapack_int( M    );                         \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _K     = detail::to_scalapack_int( K    );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _LWORK = detail::to_scalapack_int( LWORK );                        \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( &_M, &_N, &_K, A, &_IA, &_JA, _DESCA.data(), TAU, WORK, &_LWORK, \
         &INFO );                                                         \
  return INFO;                                                            \
                                                                          \
}

porgqr_impl( float,    psorgqr_ );
porgqr_impl( double,   pdorgqr_ );


}
}

