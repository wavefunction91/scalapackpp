/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/factorizations/potrf.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

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

namespace scalapackpp::wrappers {

template <>
scalapack_int
  ppotrf( const char* UPLO, scalapack_int N, float* A, scalapack_int IA, 
          scalapack_int JA, const scalapack_desc& DESCA ){


  scalapack_int INFO;
  pspotrf_( UPLO, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ppotrf( const char* UPLO, scalapack_int N, double* A, scalapack_int IA, 
          scalapack_int JA, const scalapack_desc& DESCA ){


  scalapack_int INFO;
  pdpotrf_( UPLO, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ppotrf( const char* UPLO, scalapack_int N, scomplex* A, scalapack_int IA, 
          scalapack_int JA, const scalapack_desc& DESCA ){


  scalapack_int INFO;
  pcpotrf_( UPLO, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ppotrf( const char* UPLO, scalapack_int N, dcomplex* A, scalapack_int IA, 
          scalapack_int JA, const scalapack_desc& DESCA ){


  scalapack_int INFO;
  pzpotrf_( UPLO, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}

}
