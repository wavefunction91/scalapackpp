/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/linear_systems/trtrs.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void pstrtrs_( const char* UPLO, const char* TRANS, const char* DIAG,
    const scalapack_int* N, const scalapack_int* NRHS, 
    const float* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
          float* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pdtrtrs_( const char* UPLO, const char* TRANS, const char* DIAG,
    const scalapack_int* N, const scalapack_int* NRHS, 
    const double* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
          double* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pctrtrs_( const char* UPLO, const char* TRANS, const char* DIAG,
    const scalapack_int* N, const scalapack_int* NRHS, 
    const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
          scomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pztrtrs_( const char* UPLO, const char* TRANS, const char* DIAG,
    const scalapack_int* N, const scalapack_int* NRHS, 
    const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
          dcomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );

}

namespace scalapackpp::wrappers {

template <>
scalapack_int
  ptrtrs( const char* UPLO, const char* TRANS, const char* DIAG,
    scalapack_int N, scalapack_int NRHS, 
    const float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          float* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  scalapack_int INFO;
  pstrtrs_( UPLO, TRANS, DIAG, &N, &NRHS, A, &IA, &JA, DESCA.data(),
            B, &IB, &JB, DESCB.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ptrtrs( const char* UPLO, const char* TRANS, const char* DIAG,
    scalapack_int N, scalapack_int NRHS, 
    const double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          double* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  scalapack_int INFO;
  pdtrtrs_( UPLO, TRANS, DIAG, &N, &NRHS, A, &IA, &JA, DESCA.data(),
            B, &IB, &JB, DESCB.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ptrtrs( const char* UPLO, const char* TRANS, const char* DIAG,
    scalapack_int N, scalapack_int NRHS, 
    const scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          scomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  scalapack_int INFO;
  pctrtrs_( UPLO, TRANS, DIAG, &N, &NRHS, A, &IA, &JA, DESCA.data(),
            B, &IB, &JB, DESCB.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ptrtrs( const char* UPLO, const char* TRANS, const char* DIAG,
    scalapack_int N, scalapack_int NRHS, 
    const dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          dcomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  scalapack_int INFO;
  pztrtrs_( UPLO, TRANS, DIAG, &N, &NRHS, A, &IA, &JA, DESCA.data(),
            B, &IB, &JB, DESCB.data(), &INFO );
  return INFO;

}
          
}
