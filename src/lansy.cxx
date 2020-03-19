#include <scalapackpp/wrappers/matrix_norm/lansy.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

float pslansy_( const char* NORM, const char* UPLO, const scalapack_int* N,
                const float* A, const scalapack_int* IA, const scalapack_int* JA,
                const scalapack_int* DESCA, float* WORK );
double pdlansy_( const char* NORM, const char* UPLO, const scalapack_int* N,
                 const double* A, const scalapack_int* IA, const scalapack_int* JA,
                 const scalapack_int* DESCA, double* WORK );
float pclansy_( const char* NORM, const char* UPLO, const scalapack_int* N,
                const scomplex* A, const scalapack_int* IA, const scalapack_int* JA,
                const scalapack_int* DESCA, float* WORK );
double pzlansy_( const char* NORM, const char* UPLO, const scalapack_int* N,
                 const dcomplex* A, const scalapack_int* IA, const scalapack_int* JA,
                 const scalapack_int* DESCA, double* WORK );

}

namespace scalapackpp::wrappers {

template <>
float plansy( const char* NORM, const char* UPLO, scalapack_int N, const float* A, 
              scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
              float* WORK ) {

  return pslansy_( NORM, UPLO, &N, A, &IA, &JA, DESCA.data(), WORK );

}

template <>
double plansy( const char* NORM, const char* UPLO, scalapack_int N, const double* A, 
               scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
               double* WORK ) {

  return pdlansy_( NORM, UPLO, &N, A, &IA, &JA, DESCA.data(), WORK );

}

template <>
float plansy( const char* NORM, const char* UPLO, scalapack_int N, const scomplex* A, 
              scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
              float* WORK ) {

  return pclansy_( NORM, UPLO, &N, A, &IA, &JA, DESCA.data(), WORK );

}

template <>
double plansy( const char* NORM, const char* UPLO, scalapack_int N, const dcomplex* A, 
               scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
               double* WORK ) {

  return pzlansy_( NORM, UPLO, &N, A, &IA, &JA, DESCA.data(), WORK );

}

}
