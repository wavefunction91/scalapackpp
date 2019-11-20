#include <scalapackpp/wrappers/linear_systems/getrs.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void psgetrs_( const char* TRANS, const scalapack_int* N, const scalapack_int* NRHS, 
    const float* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    float* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pdgetrs_( const char* TRANS, const scalapack_int* N, const scalapack_int* NRHS, 
    const double* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    double* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pcgetrs_( const char* TRANS, const scalapack_int* N, const scalapack_int* NRHS, 
    const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    scomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );
void pzgetrs_( const char* TRANS, const scalapack_int* N, const scalapack_int* NRHS, 
    const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
    const scalapack_int* IPIV,
    dcomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB,
    scalapack_int* INFO );

}


namespace scalapackpp::wrappers {

template <>
scalapack_int
  pgetrs( const char* TRANS, scalapack_int N, scalapack_int NRHS, 
    const float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    const scalapack_int* IPIV,
    float* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

    scalapack_int INFO;
    psgetrs_( TRANS, &N, &NRHS, A, &IA, &JA, DESCA.data(), IPIV, B, &IB, &JB, 
      DESCB.data(), &INFO );
    return INFO;
}
template <>
scalapack_int
  pgetrs( const char* TRANS, scalapack_int N, scalapack_int NRHS, 
    const double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    const scalapack_int* IPIV,
    double* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

    scalapack_int INFO;
    pdgetrs_( TRANS, &N, &NRHS, A, &IA, &JA, DESCA.data(), IPIV, B, &IB, &JB, 
      DESCB.data(), &INFO );
    return INFO;
}
template <>
scalapack_int
  pgetrs( const char* TRANS, scalapack_int N, scalapack_int NRHS, 
    const scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    const scalapack_int* IPIV,
    scomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

    scalapack_int INFO;
    pcgetrs_( TRANS, &N, &NRHS, A, &IA, &JA, DESCA.data(), IPIV, B, &IB, &JB, 
      DESCB.data(), &INFO );
    return INFO;
}
template <>
scalapack_int
  pgetrs( const char* TRANS, scalapack_int N, scalapack_int NRHS, 
    const dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    const scalapack_int* IPIV,
    dcomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

    scalapack_int INFO;
    pzgetrs_( TRANS, &N, &NRHS, A, &IA, &JA, DESCA.data(), IPIV, B, &IB, &JB, 
      DESCB.data(), &INFO );
    return INFO;
}
          
}
