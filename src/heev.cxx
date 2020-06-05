/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/eigenvalue_problem/heev.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void pcheev_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, float* W,
         scomplex* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, scomplex* WORK, const scalapack_int *LWORK,
         float* RWORK, const scalapack_int *LRWORK, scalapack_int* INFO );

void pzheev_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, double* W,
         dcomplex* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, dcomplex* WORK, const scalapack_int *LWORK,
         double* RWORK, const scalapack_int *LRWORK, scalapack_int* INFO );

void pcheevd_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, float* W,
         scomplex* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, scomplex* WORK, const scalapack_int *LWORK,
         float* RWORK, const scalapack_int* LRWORK,
         scalapack_int* IWORK, const scalapack_int* LIWORK, scalapack_int* INFO );

void pzheevd_( const char* JOBZ, const char* UPLO, const scalapack_int* N,
         dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA, double* W,
         dcomplex* Z, const scalapack_int* IZ, const scalapack_int* JZ, 
         const scalapack_int* DESCZ, dcomplex* WORK, const scalapack_int *LWORK,
         double* RWORK, const scalapack_int* LRWORK,
         scalapack_int* IWORK, const scalapack_int* LIWORK, scalapack_int* INFO );

}


namespace scalapackpp::wrappers {

template <>
scalapack_int
  pheev( const char* JOBZ, const char* UPLO, scalapack_int N,
         dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         double* W,
         dcomplex* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
         dcomplex* WORK, scalapack_int LWORK, double* RWORK, scalapack_int LRWORK ) {

  scalapack_int INFO;
  pzheev_( JOBZ, UPLO, &N, A, &IA, &JA, DESCA.data(), W, Z, &IZ, &JZ, DESCZ.data(),
           WORK, &LWORK, RWORK, &LRWORK, &INFO );
  return INFO;

}

template <>
scalapack_int
  pheev( const char* JOBZ, const char* UPLO, scalapack_int N,
         scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         float* W,
         scomplex* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
         scomplex* WORK, scalapack_int LWORK, float* RWORK, scalapack_int LRWORK ) {

  scalapack_int INFO;
  pcheev_( JOBZ, UPLO, &N, A, &IA, &JA, DESCA.data(), W, Z, &IZ, &JZ, DESCZ.data(),
           WORK, &LWORK, RWORK, &LRWORK, &INFO );
  return INFO;

}











template <>
scalapack_int
  pheevd( const char* JOBZ, const char* UPLO, scalapack_int N,
          dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          double* W,
          dcomplex* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
          dcomplex* WORK, scalapack_int LWORK, double* RWORK, scalapack_int LRWORK,
          scalapack_int* IWORK, scalapack_int LIWORK ) {

  scalapack_int INFO;
  pzheevd_( JOBZ, UPLO, &N, A, &IA, &JA, DESCA.data(), W, Z, &IZ, &JZ, DESCZ.data(),
            WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO );
  return INFO;

}

template <>
scalapack_int
  pheevd( const char* JOBZ, const char* UPLO, scalapack_int N,
          scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          float* W,
          scomplex* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
          scomplex* WORK, scalapack_int LWORK, float* RWORK, scalapack_int LRWORK,
          scalapack_int* IWORK, scalapack_int LIWORK ) {

  scalapack_int INFO;
  pcheevd_( JOBZ, UPLO, &N, A, &IA, &JA, DESCA.data(), W, Z, &IZ, &JZ, DESCZ.data(),
            WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO );
  return INFO;

}

}

