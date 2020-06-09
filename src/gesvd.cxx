/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/gesvd.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

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

namespace scalapackpp::wrappers {

template <>
scalapack_int
  pgesvd( const char* JOBU, const char* JOBVT, scalapack_int M, scalapack_int N,
         float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         float* S, float* U,  scalapack_int IU,  scalapack_int JU,  
	 const scalapack_desc& DESCU, float* VT, scalapack_int IVT, 
	 scalapack_int JVT, const scalapack_desc& DESCVT, float* WORK, 
	 scalapack_int LWORK) {

  scalapack_int INFO = 0;
  psgesvd_(JOBU, JOBVT, &M, &N, A, &IA, &JA, DESCA.data(), S, U, &IU, &JU,
    DESCU.data(), VT, &IVT, &JVT, DESCVT.data(), WORK, &LWORK, &INFO );
  return INFO;

} 

template <>
scalapack_int
  pgesvd( const char* JOBU, const char* JOBVT, scalapack_int M, scalapack_int N,
         double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         double* S, double* U,  scalapack_int IU,  scalapack_int JU,  
	 const scalapack_desc& DESCU, double* VT, scalapack_int IVT, 
	 scalapack_int JVT, const scalapack_desc& DESCVT, double* WORK, 
	 scalapack_int LWORK) {

  scalapack_int INFO = 0;
  pdgesvd_(JOBU, JOBVT, &M, &N, A, &IA, &JA, DESCA.data(), S, U, &IU, &JU,
    DESCU.data(), VT, &IVT, &JVT, DESCVT.data(), WORK, &LWORK, &INFO );
  return INFO;

} 

template <>
scalapack_int
  pgesvd( const char* JOBU, const char* JOBVT, scalapack_int M, scalapack_int N,
         scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         float* S, scomplex* U,  scalapack_int IU,  scalapack_int JU,  
	 const scalapack_desc& DESCU, scomplex* VT, scalapack_int IVT, 
	 scalapack_int JVT, const scalapack_desc& DESCVT, scomplex* WORK, 
	 scalapack_int LWORK, float* RWORK) {

  scalapack_int INFO=0;
  pcgesvd_(JOBU, JOBVT, &M, &N, A, &IA, &JA, DESCA.data(), S, U, &IU, &JU,
    DESCU.data(), VT, &IVT, &JVT, DESCVT.data(), WORK, &LWORK, RWORK, &INFO );
  return INFO;

} 


template <>
scalapack_int
  pgesvd( const char* JOBU, const char* JOBVT, scalapack_int M, scalapack_int N,
         dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         double* S, dcomplex* U,  scalapack_int IU,  scalapack_int JU,  
	 const scalapack_desc& DESCU, dcomplex* VT, scalapack_int IVT, 
	 scalapack_int JVT, const scalapack_desc& DESCVT, dcomplex* WORK, 
	 scalapack_int LWORK, double* RWORK) {

  scalapack_int INFO=0;
  pzgesvd_(JOBU, JOBVT, &M, &N, A, &IA, &JA, DESCA.data(), S, U, &IU, &JU,
    DESCU.data(), VT, &IVT, &JVT, DESCVT.data(), WORK, &LWORK, RWORK, &INFO );
  return INFO;

} 


}
