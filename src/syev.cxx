/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/eigenvalue_problem/syev.hpp>

using scalapackpp::scalapack_int;

// Prototypes
extern "C" {

void pssyev_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         float* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, float* W,
         float* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, float* WORK, const scalapack_int *LWORK,
         scalapack_int* INFO );

void pdsyev_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         double* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, double* W,
         double* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, double* WORK, const scalapack_int *LWORK,
         scalapack_int* INFO );

void pssyevd_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         float* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, float* W,
         float* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, float* WORK, const scalapack_int *LWORK,
         scalapack_int* IWORK, const scalapack_int* LIWORK, scalapack_int* INFO );

void pdsyevd_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         double* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, double* W,
         double* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, double* WORK, const scalapack_int *LWORK,
         scalapack_int* IWORK, const scalapack_int* LIWORK, scalapack_int* INFO );

}


namespace scalapackpp::wrappers {

template <>
scalapack_int
  psyev( const char* JOBZ, const char* UPLO, scalapack_int N,
         double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         double* W,
         double* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
         double* WORK, scalapack_int LWORK ) {

  scalapack_int INFO;
  pdsyev_( JOBZ, UPLO, &N, A, &IA, &JA, DESCA.data(), W, Z, &IZ, &JZ, DESCZ.data(),
           WORK, &LWORK, &INFO );
  return INFO;

}

template <>
scalapack_int
  psyev( const char* JOBZ, const char* UPLO, scalapack_int N,
         float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         float* W,
         float* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
         float* WORK, scalapack_int LWORK ) {

  scalapack_int INFO;
  pssyev_( JOBZ, UPLO, &N, A, &IA, &JA, DESCA.data(), W, Z, &IZ, &JZ, DESCZ.data(),
           WORK, &LWORK, &INFO );
  return INFO;

}











template <>
scalapack_int
  psyevd( const char* JOBZ, const char* UPLO, scalapack_int N,
          double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          double* W,
          double* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
          double* WORK, scalapack_int LWORK, scalapack_int* IWORK, 
          scalapack_int LIWORK ) {

  scalapack_int INFO;
  pdsyevd_( JOBZ, UPLO, &N, A, &IA, &JA, DESCA.data(), W, Z, &IZ, &JZ, DESCZ.data(),
            WORK, &LWORK, IWORK, &LIWORK, &INFO );
  return INFO;

}

template <>
scalapack_int
  psyevd( const char* JOBZ, const char* UPLO, scalapack_int N,
          float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          float* W,
          float* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
          float* WORK, scalapack_int LWORK, scalapack_int* IWORK, 
          scalapack_int LIWORK ) {

  scalapack_int INFO;
  pssyevd_( JOBZ, UPLO, &N, A, &IA, &JA, DESCA.data(), W, Z, &IZ, &JZ, DESCZ.data(),
            WORK, &LWORK, IWORK, &LIWORK, &INFO );
  return INFO;

}

}

