#include <scalapackpp/wrappers/matrix_norm/lange.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

float pslange_( const char* NORM, const scalapack_int* M, const scalapack_int* N,
                const float* A, const scalapack_int* IA, const scalapack_int* JA,
                const scalapack_int* DESCA, float* WORK );
double pdlange_( const char* NORM, const scalapack_int* M, const scalapack_int* N,
                 const double* A, const scalapack_int* IA, const scalapack_int* JA,
                 const scalapack_int* DESCA, double* WORK );
float pclange_( const char* NORM, const scalapack_int* M, const scalapack_int* N,
                const scomplex* A, const scalapack_int* IA, const scalapack_int* JA,
                const scalapack_int* DESCA, float* WORK );
double pzlange_( const char* NORM, const scalapack_int* M, const scalapack_int* N,
                 const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA,
                 const scalapack_int* DESCA, double* WORK );

}

namespace scalapackpp::wrappers {

template <>
float plange( const char* NORM, scalapack_int M, scalapack_int N, const float* A, 
              scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
              float* WORK ) {

  return pslange_( NORM, &M, &N, A, &IA, &JA, DESCA.data(), WORK );

}

template <>
double plange( const char* NORM, scalapack_int M, scalapack_int N, const double* A, 
               scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
               double* WORK ) {

  return pdlange_( NORM, &M, &N, A, &IA, &JA, DESCA.data(), WORK );

}

template <>
float plange( const char* NORM, scalapack_int M, scalapack_int N, const scomplex* A, 
              scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
              float* WORK ) {

  return pclange_( NORM, &M, &N, A, &IA, &JA, DESCA.data(), WORK );

}

template <>
double plange( const char* NORM, scalapack_int M, scalapack_int N, const dcomplex* A, 
               scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
               double* WORK ) {

  return pzlange_( NORM, &M, &N, A, &IA, &JA, DESCA.data(), WORK );

}

}
