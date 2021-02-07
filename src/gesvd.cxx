/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/gesvd.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;
using scalapackpp::internal::dcomplex;
using scalapackpp::internal::scomplex;

// Prototypes
extern "C" {

void psgesvd_(const char* JOBU, const char* JOBVT, const scalapack_int* M, 
	      const scalapack_int* N, float* A, const scalapack_int* IA, 
	      const scalapack_int* JA, const scalapack_int* DESCA, float* S,
              float* U, const scalapack_int* IU, const scalapack_int* JU,  
	      const scalapack_int* DESCU, float* VT, const scalapack_int* IVT, 
	      const scalapack_int* JVT, const scalapack_int* DESCVT,
	      float* WORK, const scalapack_int* LWORK, scalapack_int* INFO); 

void pdgesvd_(const char* JOBU, const char* JOBVT, const scalapack_int* M, 
	      const scalapack_int* N, double* A, const scalapack_int* IA, 
	      const scalapack_int* JA, const scalapack_int* DESCA, double* S,
              double* U, const scalapack_int* IU, const scalapack_int* JU,  
	      const scalapack_int* DESCU, double* VT, const scalapack_int* IVT, 
	      const scalapack_int* JVT, const scalapack_int* DESCVT,
	      double* WORK, const scalapack_int* LWORK, scalapack_int* INFO); 

void pcgesvd_(const char* JOBU, const char* JOBVT, const scalapack_int* M, 
	      const scalapack_int* N, scomplex* A, const scalapack_int* IA, 
	      const scalapack_int* JA, const scalapack_int* DESCA, float* S,
              scomplex* U, const scalapack_int* IU, const scalapack_int* JU,  
	      const scalapack_int* DESCU, scomplex* VT, const scalapack_int* IVT, 
	      const scalapack_int* JVT, const scalapack_int* DESCVT,
	      scomplex* WORK, const scalapack_int* LWORK, float* RWORK, 
	      scalapack_int* INFO); 

void pzgesvd_(const char* JOBU, const char* JOBVT, const scalapack_int* M, 
	      const scalapack_int* N, dcomplex* A, const scalapack_int* IA, 
	      const scalapack_int* JA, const scalapack_int* DESCA, double* S,
              dcomplex* U, const scalapack_int* IU, const scalapack_int* JU,  
	      const scalapack_int* DESCU, dcomplex* VT, const scalapack_int* IVT, 
	      const scalapack_int* JVT, const scalapack_int* DESCVT,
	      dcomplex* WORK, const scalapack_int* LWORK, double* RWORK,
	      scalapack_int* INFO); 

}

namespace scalapackpp {
namespace wrappers    {

#define pgesvd_real_impl(F,FNAME)                                         \
template <>                                                               \
int64_t                                                                   \
  pgesvd( const char* JOBU, const char* JOBVT, int64_t M, int64_t N,      \
         F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,       \
         detail::real_t<F>* S, F* U,  int64_t IU,  int64_t JU,            \
	 const scalapack_desc& DESCU, F* VT, int64_t IVT,                       \
	 int64_t JVT, const scalapack_desc& DESCVT, F* WORK,                    \
	 int64_t LWORK) {                                                       \
                                                                          \
  auto _M     = detail::to_scalapack_int( M    );                         \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
  auto _IU    = detail::to_scalapack_int( IU   );                         \
  auto _JU    = detail::to_scalapack_int( JU   );                         \
  auto _DESCU = detail::to_scalapack_int( DESCU );                        \
  auto _IVT    = detail::to_scalapack_int( IVT   );                       \
  auto _JVT    = detail::to_scalapack_int( JVT   );                       \
  auto _DESCVT = detail::to_scalapack_int( DESCVT );                      \
  auto _LWORK  = detail::to_scalapack_int( LWORK  );                      \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME(JOBU, JOBVT, &_M, &_N, A, &_IA, &_JA, _DESCA.data(),              \
        S, U, &_IU, &_JU, _DESCU.data(), VT, &_IVT, &_JVT, _DESCVT.data(),\
        WORK, &_LWORK, &INFO );                                           \
  return INFO;                                                            \
                                                                          \
} 

pgesvd_real_impl( float,  psgesvd_ );
pgesvd_real_impl( double, pdgesvd_ );

#define pgesvd_complex_impl(F,FNAME)                                      \
template <>                                                               \
int64_t                                                                   \
  pgesvd( const char* JOBU, const char* JOBVT, int64_t M, int64_t N,      \
         F* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,       \
         detail::real_t<F>* S, F* U,  int64_t IU,  int64_t JU,            \
	 const scalapack_desc& DESCU, F* VT, int64_t IVT,                       \
	 int64_t JVT, const scalapack_desc& DESCVT, F* WORK,                    \
	 int64_t LWORK, detail::real_t<F>* RWORK) {                             \
                                                                          \
  auto _M     = detail::to_scalapack_int( M    );                         \
  auto _N     = detail::to_scalapack_int( N    );                         \
  auto _IA    = detail::to_scalapack_int( IA   );                         \
  auto _JA    = detail::to_scalapack_int( JA   );                         \
  auto _DESCA = detail::to_scalapack_int( DESCA );                        \
  auto _IU    = detail::to_scalapack_int( IU   );                         \
  auto _JU    = detail::to_scalapack_int( JU   );                         \
  auto _DESCU = detail::to_scalapack_int( DESCU );                        \
  auto _IVT    = detail::to_scalapack_int( IVT   );                       \
  auto _JVT    = detail::to_scalapack_int( JVT   );                       \
  auto _DESCVT = detail::to_scalapack_int( DESCVT );                      \
  auto _LWORK  = detail::to_scalapack_int( LWORK  );                      \
                                                                          \
  scalapack_int INFO;                                                     \
  FNAME(JOBU, JOBVT, &_M, &_N, A, &_IA, &_JA, _DESCA.data(),              \
        S, U, &_IU, &_JU, _DESCU.data(), VT, &_IVT, &_JVT, _DESCVT.data(),\
        WORK, &_LWORK, RWORK, &INFO );                                    \
  return INFO;                                                            \
                                                                          \
} 

pgesvd_complex_impl( scomplex, pcgesvd_ );
pgesvd_complex_impl( dcomplex, pzgesvd_ );

}
}
