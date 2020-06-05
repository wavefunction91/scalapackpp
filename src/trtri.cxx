/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/matrix_inverse/trtri.hpp>
#include <iostream>


using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

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

template <>
scalapack_int
  ptrtri( const char* UPLO, const char* DIAG, scalapack_int N,
          float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA ) {


  scalapack_int INFO;
  pstrtri_( UPLO, DIAG, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ptrtri( const char* UPLO, const char* DIAG, scalapack_int N,
          double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA ) {

  scalapack_int INFO;
  pdtrtri_( UPLO, DIAG, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ptrtri( const char* UPLO, const char* DIAG, scalapack_int N,
          scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA ) {

  scalapack_int INFO;
  pctrtri_( UPLO, DIAG, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ptrtri( const char* UPLO, const char* DIAG, scalapack_int N,
          dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA ) {

  std::cout << *UPLO << ", " << *DIAG << ", " << N << ", " << A << ", "
            << IA << ", " << JA << std::endl;
  scalapack_int INFO;
  pztrtri_( UPLO, DIAG, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}


}
