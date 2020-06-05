/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/geadd.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void psgeadd_( const char* TRANS, const scalapack_int* M, const scalapack_int* N, 
               const float* ALPHA, const float* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA,
               const float* BETA, float* C, const scalapack_int* IC, 
               const scalapack_int* JC, const scalapack_int* DESCC );
void pdgeadd_( const char* TRANS, const scalapack_int* M, const scalapack_int* N, 
               const double* ALPHA, const double* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA,
               const double* BETA, double* C, const scalapack_int* IC, 
               const scalapack_int* JC, const scalapack_int* DESCC );
void pcgeadd_( const char* TRANS, const scalapack_int* M, const scalapack_int* N, 
               const scomplex* ALPHA, const scomplex* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA,
               const scomplex* BETA, scomplex* C, const scalapack_int* IC, 
               const scalapack_int* JC, const scalapack_int* DESCC );
void pzgeadd_( const char* TRANS, const scalapack_int* M, const scalapack_int* N, 
               const dcomplex* ALPHA, const dcomplex* A, const scalapack_int* IA, 
               const scalapack_int* JA, const scalapack_int* DESCA,
               const dcomplex* BETA, dcomplex* C, const scalapack_int* IC, 
               const scalapack_int* JC, const scalapack_int* DESCC );

}



namespace scalapackpp::wrappers {

template <>
void pgeadd( const char* TRANS, scalapack_int M, scalapack_int N, float ALPHA,
             const float* A, scalapack_int IA, scalapack_int JA, 
             const scalapack_desc& DESCA,
             float BETA, 
             float* C, scalapack_int IC, scalapack_int JC, 
             const scalapack_desc& DESCC ) {


  psgeadd_( TRANS, &M, &N, &ALPHA, A, &IA, &JA, DESCA.data(), &BETA, C, &IC, &JC,
            DESCC.data() );

}

template <>
void pgeadd( const char* TRANS, scalapack_int M, scalapack_int N, double ALPHA,
             const double* A, scalapack_int IA, scalapack_int JA, 
             const scalapack_desc& DESCA,
             double BETA, 
             double* C, scalapack_int IC, scalapack_int JC, 
             const scalapack_desc& DESCC ) {


  pdgeadd_( TRANS, &M, &N, &ALPHA, A, &IA, &JA, DESCA.data(), &BETA, C, &IC, &JC,
            DESCC.data() );

}

template <>
void pgeadd( const char* TRANS, scalapack_int M, scalapack_int N, scomplex ALPHA,
             const scomplex* A, scalapack_int IA, scalapack_int JA, 
             const scalapack_desc& DESCA,
             scomplex BETA, 
             scomplex* C, scalapack_int IC, scalapack_int JC, 
             const scalapack_desc& DESCC ) {


  pcgeadd_( TRANS, &M, &N, &ALPHA, A, &IA, &JA, DESCA.data(), &BETA, C, &IC, &JC,
            DESCC.data() );

}

template <>
void pgeadd( const char* TRANS, scalapack_int M, scalapack_int N, dcomplex ALPHA,
             const dcomplex* A, scalapack_int IA, scalapack_int JA, 
             const scalapack_desc& DESCA,
             dcomplex BETA, 
             dcomplex* C, scalapack_int IC, scalapack_int JC, 
             const scalapack_desc& DESCC ) {


  pzgeadd_( TRANS, &M, &N, &ALPHA, A, &IA, &JA, DESCA.data(), &BETA, C, &IC, &JC,
            DESCC.data() );

}

}
