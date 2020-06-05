/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/linear_systems/gesv.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void psgesv_( const scalapack_int* N, const scalapack_int* NRHS, 
    float* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    float* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pdgesv_( const scalapack_int* N, const scalapack_int* NRHS, 
    double* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    double* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pcgesv_( const scalapack_int* N, const scalapack_int* NRHS, 
    scomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    scomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pzgesv_( const scalapack_int* N, const scalapack_int* NRHS, 
    dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    dcomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );

}


namespace scalapackpp::wrappers {

template <>
scalapack_int
  pgesv( scalapack_int N, scalapack_int NRHS, 
    float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    const scalapack_int* IPIV,
    float* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

    scalapack_int INFO;
    psgesv_( &N, &NRHS, A, &IA, &JA, DESCA.data(), IPIV, B, &IB, &JB, 
      DESCB.data(), &INFO );
    return INFO;
}
template <>
scalapack_int
  pgesv( scalapack_int N, scalapack_int NRHS, 
    double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    const scalapack_int* IPIV,
    double* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

    scalapack_int INFO;
    pdgesv_( &N, &NRHS, A, &IA, &JA, DESCA.data(), IPIV, B, &IB, &JB, 
      DESCB.data(), &INFO );
    return INFO;
}
template <>
scalapack_int
  pgesv( scalapack_int N, scalapack_int NRHS, 
    scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    const scalapack_int* IPIV,
    scomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

    scalapack_int INFO;
    pcgesv_( &N, &NRHS, A, &IA, &JA, DESCA.data(), IPIV, B, &IB, &JB, 
      DESCB.data(), &INFO );
    return INFO;
}
template <>
scalapack_int
  pgesv( scalapack_int N, scalapack_int NRHS, 
    dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    const scalapack_int* IPIV,
    dcomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

    scalapack_int INFO;
    pzgesv_( &N, &NRHS, A, &IA, &JA, DESCA.data(), IPIV, B, &IB, &JB, 
      DESCB.data(), &INFO );
    return INFO;
}
          
}
