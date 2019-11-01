#include <scalapackpp/wrappers/gemr2d.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;


extern "C" {


void psgemr2d_( const scalapack_int* M, const scalapack_int* N, 
    const float* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    float * B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, const scalapack_int* ICONTEXT );
void pdgemr2d_( const scalapack_int* M, const scalapack_int* N, 
    const double* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    double * B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, const scalapack_int* ICONTEXT );
void pcgemr2d_( const scalapack_int* M, const scalapack_int* N, 
    const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    scomplex * B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, const scalapack_int* ICONTEXT );
void pzgemr2d_( const scalapack_int* M, const scalapack_int* N, 
    const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    dcomplex * B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, const scalapack_int* ICONTEXT );

}


namespace scalapackpp::wrappers {

template <>
void pgemr2d( scalapack_int M, scalapack_int N, 
    const float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    float* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
    scalapack_int ICONTEXT ) {

  psgemr2d_( &M, &N, A, &IA, &JA, DESCA.data(), B, &IB, &JB, DESCB.data(), 
             &ICONTEXT );

}



template <>
void pgemr2d( scalapack_int M, scalapack_int N, 
    const double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    double* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
    scalapack_int ICONTEXT ) {

  pdgemr2d_( &M, &N, A, &IA, &JA, DESCA.data(), B, &IB, &JB, DESCB.data(), 
             &ICONTEXT );

}



template <>
void pgemr2d( scalapack_int M, scalapack_int N, 
    const scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    scomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
    scalapack_int ICONTEXT ) {

  pcgemr2d_( &M, &N, A, &IA, &JA, DESCA.data(), B, &IB, &JB, DESCB.data(), 
             &ICONTEXT );

}



template <>
void pgemr2d( scalapack_int M, scalapack_int N, 
    const dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    dcomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
    scalapack_int ICONTEXT ) {

  pzgemr2d_( &M, &N, A, &IA, &JA, DESCA.data(), B, &IB, &JB, DESCB.data(), 
             &ICONTEXT );

}

}
