/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/householder/unmqr.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pcunmqr_( const char* SIDE, const char* TRANS, const scalapack_int* M, const scalapack_int* N, 
  const scalapack_int* K, const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  const scomplex* TAU, scomplex* C, const scalapack_int* IC, const scalapack_int* JC, const scalapack_int* DESCC,
  scomplex* WORK, const scalapack_int* LWORK, scalapack_int* INFO );
void pzunmqr_( const char* SIDE, const char* TRANS, const scalapack_int* M, const scalapack_int* N, 
  const scalapack_int* K, const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  const dcomplex* TAU, dcomplex* C, const scalapack_int* IC, const scalapack_int* JC, const scalapack_int* DESCC,
  dcomplex* WORK, const scalapack_int* LWORK, scalapack_int* INFO );

}


namespace scalapackpp {
namespace wrappers    {

#define punmqr_impl(F,FNAME)                                              \
template <>                                                               \
int64_t                                                                   \
  punmqr( const char* SIDE, const char* TRANS, int64_t M, int64_t N,      \
          int64_t K, const F* A, int64_t IA, int64_t JA,                  \
          const scalapack_desc& DESCA, const F* TAU, F* C, int64_t IC,    \
          int64_t JC, const scalapack_desc& DESCC, F* WORK,               \
          int64_t LWORK ) {                                               \
                                                                          \
  auto _M     = detail::to_scalapack_int( M    );                         \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _K     = detail::to_scalapack_int( K    );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _IC    = detail::to_scalapack_int( IC   );                         \
  auto _JC    = detail::to_scalapack_int( JC   );                         \
  auto _LWORK = detail::to_scalapack_int( LWORK );                        \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
  auto _DESCC = detail::to_scalapack_int( DESCC );                        \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME( SIDE, TRANS, &_M, &_N, &_K, A, &_IA, &_JA, _DESCA.data(), TAU,   \
         C, &_IC, &_JC, _DESCC.data(), WORK, &_LWORK, &INFO );            \
  return INFO;                                                            \
                                                                          \
}

punmqr_impl( scomplex, pcunmqr_ );
punmqr_impl( dcomplex, pzunmqr_ );


}
}

