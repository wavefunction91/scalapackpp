#include <scalapackpp/wrappers/potrs.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void pspotrs_( const char* UPLO, const scalapack_int* N, const scalapack_int* NRHS, 
    const float* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    float* B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, scalapack_int* INFO );
void pdpotrs_( const char* UPLO, const scalapack_int* N, const scalapack_int* NRHS, 
    const double* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    double* B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, scalapack_int* INFO );
void pcpotrs_( const char* UPLO, const scalapack_int* N, const scalapack_int* NRHS, 
    const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    scomplex* B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, scalapack_int* INFO );
void pzpotrs_( const char* UPLO, const scalapack_int* N, const scalapack_int* NRHS, 
    const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
    const scalapack_int* DESCA,
    dcomplex* B, const scalapack_int* IB, const scalapack_int* JB, 
    const scalapack_int* DESCB, scalapack_int* INFO );

}




namespace scalapackpp::wrappers {

template <>
scalapack_int
  ppotrs( const char* UPLO, scalapack_int N, scalapack_int NRHS, 
    const float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    float* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  scalapack_int INFO;
  pspotrs_( UPLO, &N, &NRHS, A, &IA, &JA, DESCA.data(), B, &IB, &JB, DESCB.data(),
            &INFO );
  return INFO;

}
template <>
scalapack_int
  ppotrs( const char* UPLO, scalapack_int N, scalapack_int NRHS, 
    const double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    double* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  scalapack_int INFO;
  pdpotrs_( UPLO, &N, &NRHS, A, &IA, &JA, DESCA.data(), B, &IB, &JB, DESCB.data(),
            &INFO );
  return INFO;

}

template <>
scalapack_int
  ppotrs( const char* UPLO, scalapack_int N, scalapack_int NRHS, 
    const scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    scomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  scalapack_int INFO;
  pcpotrs_( UPLO, &N, &NRHS, A, &IA, &JA, DESCA.data(), B, &IB, &JB, DESCB.data(),
            &INFO );
  return INFO;

}

template <>
scalapack_int
  ppotrs( const char* UPLO, scalapack_int N, scalapack_int NRHS, 
    const dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    dcomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  scalapack_int INFO;
  pzpotrs_( UPLO, &N, &NRHS, A, &IA, &JA, DESCA.data(), B, &IB, &JB, DESCB.data(),
            &INFO );
  return INFO;

}

}
