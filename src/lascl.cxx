/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/lascl.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void pslascl_( const char*, const float*, const float*, const scalapack_int* M, 
               const scalapack_int* N, float* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA ); 

void pdlascl_( const char*, const double*, const double*, const scalapack_int* M, 
               const scalapack_int* N, double* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA ); 

void pclascl_( const char*, const float*, const float*, const scalapack_int* M, 
               const scalapack_int* N, scomplex* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA ); 

void pzlascl_( const char*, const double*, const double*, const scalapack_int* M, 
               const scalapack_int* N, dcomplex* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA ); 

}

namespace scalapackpp {
namespace wrappers    {

#define plascl_impl(F, FNAME)                                                  \
template <>                                                                    \
void plascl( const char* TYPE, detail::real_t<F> CTO, detail::real_t<F> CFROM, \
             int64_t M, int64_t N, F* A, int64_t IA, int64_t JA,               \
             const scalapack_desc& DESCA) {                                    \
                                                                               \
  auto _M     = detail::to_scalapack_int( M  );                                \
  auto _N     = detail::to_scalapack_int( N  );                                \
  auto _IA    = detail::to_scalapack_int( IA );                                \
  auto _JA    = detail::to_scalapack_int( JA );                                \
                                                                               \
  auto _DESCA = detail::to_scalapack_int( DESCA );                             \
                                                                               \
  FNAME( TYPE, &CTO, &CFROM, &_M, &_N, A, &_IA, &_JA, _DESCA.data());          \
}

plascl_impl( float,    pslascl_ );
plascl_impl( double,   pdlascl_ );
plascl_impl( scomplex, pclascl_ );
plascl_impl( dcomplex, pzlascl_ );

}
}

