#include <scalapackpp/wrappers/pblas/gemm.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void psgemm_( const char* TRANSA, const char* TRANSB,
         const scalapack_int* M, const scalapack_int* N, const scalapack_int* K, 
         const float* ALPHA, 
         const float* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const float* B, const scalapack_int* IB, const scalapack_int* JB, 
         const scalapack_int* DESCB,
         const float* BETA,
         float* C, const scalapack_int* IC, const scalapack_int* JC, 
         const scalapack_int* DESCC );
void pdgemm_( const char* TRANSA, const char* TRANSB,
         const scalapack_int* M, const scalapack_int* N, const scalapack_int* K, 
         const double* ALPHA, 
         const double* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const double* B, const scalapack_int* IB, const scalapack_int* JB, 
         const scalapack_int* DESCB,
         const double* BETA,
         double* C, const scalapack_int* IC, const scalapack_int* JC, 
         const scalapack_int* DESCC );
void pcgemm_( const char* TRANSA, const char* TRANSB,
         const scalapack_int* M, const scalapack_int* N, const scalapack_int* K, 
         const scomplex* ALPHA, 
         const scomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const scomplex* B, const scalapack_int* IB, const scalapack_int* JB, 
         const scalapack_int* DESCB,
         const scomplex* BETA,
         scomplex* C, const scalapack_int* IC, const scalapack_int* JC, 
         const scalapack_int* DESCC );
void pzgemm_( const char* TRANSA, const char* TRANSB,
         const scalapack_int* M, const scalapack_int* N, const scalapack_int* K, 
         const dcomplex* ALPHA, 
         const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA, 
         const scalapack_int* DESCA,
         const dcomplex* B, const scalapack_int* IB, const scalapack_int* JB, 
         const scalapack_int* DESCB,
         const dcomplex* BETA,
         dcomplex* C, const scalapack_int* IC, const scalapack_int* JC, 
         const scalapack_int* DESCC );

}



namespace scalapackpp::wrappers {

template <>
void 
  pgemm( const char* TRANSA, const char* TRANSB,
         scalapack_int M, scalapack_int N, scalapack_int K, float ALPHA, 
         const float* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         const float* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
         float BETA,
         float* C, scalapack_int IC, scalapack_int JC, const scalapack_desc& DESCC ) {

  psgemm_( TRANSA, TRANSB, &M, &N, &K, &ALPHA, A, &IA, &JA, DESCA.data(), 
           B, &IB, &JB, DESCB.data(), &BETA, C, &IC, &JC, DESCC.data() );

}

template <>
void 
  pgemm( const char* TRANSA, const char* TRANSB,
         scalapack_int M, scalapack_int N, scalapack_int K, double ALPHA, 
         const double* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         const double* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
         double BETA,
         double* C, scalapack_int IC, scalapack_int JC, const scalapack_desc& DESCC ) {

  pdgemm_( TRANSA, TRANSB, &M, &N, &K, &ALPHA, A, &IA, &JA, DESCA.data(), 
           B, &IB, &JB, DESCB.data(), &BETA, C, &IC, &JC, DESCC.data() );

}

template <>
void 
  pgemm( const char* TRANSA, const char* TRANSB,
         scalapack_int M, scalapack_int N, scalapack_int K, scomplex ALPHA, 
         const scomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         const scomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
         scomplex BETA,
         scomplex* C, scalapack_int IC, scalapack_int JC, const scalapack_desc& DESCC ) {

  pcgemm_( TRANSA, TRANSB, &M, &N, &K, &ALPHA, A, &IA, &JA, DESCA.data(), 
           B, &IB, &JB, DESCB.data(), &BETA, C, &IC, &JC, DESCC.data() );

}



template <>
void 
  pgemm( const char* TRANSA, const char* TRANSB,
         scalapack_int M, scalapack_int N, scalapack_int K, dcomplex ALPHA, 
         const dcomplex* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         const dcomplex* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
         dcomplex BETA,
         dcomplex* C, scalapack_int IC, scalapack_int JC, const scalapack_desc& DESCC ) {

  pzgemm_( TRANSA, TRANSB, &M, &N, &K, &ALPHA, A, &IA, &JA, DESCA.data(), 
           B, &IB, &JB, DESCB.data(), &BETA, C, &IC, &JC, DESCC.data() );

}

}

