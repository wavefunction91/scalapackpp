/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/matrix_norm/lanhe.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

float pclanhe_( const char* NORM, const char* UPLO, const scalapack_int* N,
                const scomplex* A, const scalapack_int* IA, const scalapack_int* JA,
                const scalapack_int* DESCA, float* WORK );
double pzlanhe_( const char* NORM, const char* UPLO, const scalapack_int* N,
                 const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA,
                 const scalapack_int* DESCA, double* WORK );

}

namespace scalapackpp::wrappers {

template <>
float planhe( const char* NORM, const char* UPLO, scalapack_int N, const scomplex* A, 
              scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
              float* WORK ) {

  return pclanhe_( NORM, UPLO, &N, A, &IA, &JA, DESCA.data(), WORK );

}

template <>
double planhe( const char* NORM, const char* UPLO, scalapack_int N, const dcomplex* A, 
               scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
               double* WORK ) {

  return pzlanhe_( NORM, UPLO, &N, A, &IA, &JA, DESCA.data(), WORK );

}

}
