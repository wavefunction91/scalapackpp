#include <scalapackpp/wrappers/trmm.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void pstrmm_( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         const scalapack_int* M, const scalapack_int* N, const float* ALPHA, 
         const float* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
         float* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB );
void pdtrmm_( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         const scalapack_int* M, const scalapack_int* N, const double* ALPHA, 
         const double* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
         double* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB );
void pctrmm_( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         const scalapack_int* M, const scalapack_int* N, const scomplex* ALPHA, 
         const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
         scomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB );
void pztrmm_( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         const scalapack_int* M, const scalapack_int* N, const dcomplex* ALPHA, 
         const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA,
         dcomplex* B, const scalapack_int* IB, const scalapack_int* JB, const scalapack_int* DESCB );

}

namespace scalapackpp::wrappers {

template <>
void ptrmm( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         scalapack_int M, scalapack_int N, float ALPHA, 
         const float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         float* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  pstrmm_( SIDE, UPLO, TRANS, DIAG, &M, &N, &ALPHA, A, &IA, &JA, DESCA.data(),
           B, &IB, &JB, DESCB.data() );

}

template <>
void ptrmm( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         scalapack_int M, scalapack_int N, double ALPHA, 
         const double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         double* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  pdtrmm_( SIDE, UPLO, TRANS, DIAG, &M, &N, &ALPHA, A, &IA, &JA, DESCA.data(),
           B, &IB, &JB, DESCB.data() );

}

template <>
void ptrmm( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         scalapack_int M, scalapack_int N, scomplex ALPHA, 
         const scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         scomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  pctrmm_( SIDE, UPLO, TRANS, DIAG, &M, &N, &ALPHA, A, &IA, &JA, DESCA.data(),
           B, &IB, &JB, DESCB.data() );

}


template <>
void ptrmm( const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG,
         scalapack_int M, scalapack_int N, dcomplex ALPHA, 
         const dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         dcomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  pztrmm_( SIDE, UPLO, TRANS, DIAG, &M, &N, &ALPHA, A, &IA, &JA, DESCA.data(),
           B, &IB, &JB, DESCB.data() );

}

}
