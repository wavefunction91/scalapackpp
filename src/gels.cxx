/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/linear_systems/gels.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psgels_( const char* TRANS, const scalapack_int* M, const scalapack_int* N,
              const scalapack_int* NRHS, float* A, const scalapack_int* IA,
              const scalapack_int* JA, const scalapack_int* DESCA,
              float* B, const scalapack_int* IB, const scalapack_int* JB,
              const scalapack_int* DESCB, float* WORK, const scalapack_int* LWORK,
              scalapack_int* INFO );
void pdgels_( const char* TRANS, const scalapack_int* M, const scalapack_int* N,
              const scalapack_int* NRHS, double* A, const scalapack_int* IA,
              const scalapack_int* JA, const scalapack_int* DESCA,
              double* B, const scalapack_int* IB, const scalapack_int* JB,
              const scalapack_int* DESCB, double* WORK, const scalapack_int* LWORK,
              scalapack_int* INFO );
void pcgels_( const char* TRANS, const scalapack_int* M, const scalapack_int* N,
              const scalapack_int* NRHS, scomplex* A, const scalapack_int* IA,
              const scalapack_int* JA, const scalapack_int* DESCA,
              scomplex* B, const scalapack_int* IB, const scalapack_int* JB,
              const scalapack_int* DESCB, scomplex* WORK, const scalapack_int* LWORK,
              scalapack_int* INFO );
void pzgels_( const char* TRANS, const scalapack_int* M, const scalapack_int* N,
              const scalapack_int* NRHS, dcomplex* A, const scalapack_int* IA,
              const scalapack_int* JA, const scalapack_int* DESCA,
              dcomplex* B, const scalapack_int* IB, const scalapack_int* JB,
              const scalapack_int* DESCB, dcomplex* WORK, const scalapack_int* LWORK,
              scalapack_int* INFO );

}

namespace scalapackpp {
namespace wrappers    {

#define pgels_impl(F,FNAME)                                               \
template <>                                                               \
int64_t                                                                   \
  pgels( const char* TRANS, int64_t M, int64_t N, int64_t NRHS,           \
    F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,            \
    F* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB,            \
    F* WORK, int64_t LWORK ) {                                            \
                                                                          \
  auto _M     = detail::to_scalapack_int( M    );                         \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _NRHS  = detail::to_scalapack_int( NRHS );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _IB    = detail::to_scalapack_int( IB   );                         \
  auto _JB    = detail::to_scalapack_int( JB   );                         \
  auto _LWORK = detail::to_scalapack_int( LWORK );                        \
                                                                          \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
  auto _DESCB = detail::to_scalapack_int( DESCB );                        \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( TRANS, &_M, &_N, &_NRHS, A, &_IA, &_JA, _DESCA.data(),           \
         B, &_IB, &_JB, _DESCB.data(), WORK, &_LWORK, &INFO );            \
  return INFO;                                                            \
                                                                          \
}

pgels_impl( float,    psgels_ );
pgels_impl( double,   pdgels_ );
pgels_impl( scomplex, pcgels_ );
pgels_impl( dcomplex, pzgels_ );

}
}
