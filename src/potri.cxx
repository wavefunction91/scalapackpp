#include <scalapackpp/wrappers/matrix_inverse/potri.hpp>


using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void pspotri_( const char* UPLO, const scalapack_int* N, float* A, const scalapack_int* IA, 
          const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* INFO );
void pdpotri_( const char* UPLO, const scalapack_int* N, double* A, const scalapack_int* IA, 
          const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* INFO );
void pcpotri_( const char* UPLO, const scalapack_int* N, scomplex* A, const scalapack_int* IA, 
          const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* INFO );
void pzpotri_( const char* UPLO, const scalapack_int* N, dcomplex* A, const scalapack_int* IA, 
          const scalapack_int* JA, const scalapack_int* DESCA, scalapack_int* INFO );

}

namespace scalapackpp::wrappers {

template <>
scalapack_int
  ppotri( const char* UPLO, scalapack_int N, float* A, scalapack_int IA, 
          scalapack_int JA, const scalapack_desc& DESCA ) {

  scalapack_int INFO;
  pspotri_( UPLO, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ppotri( const char* UPLO, scalapack_int N, double* A, scalapack_int IA, 
          scalapack_int JA, const scalapack_desc& DESCA ) {

  scalapack_int INFO;
  pdpotri_( UPLO, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ppotri( const char* UPLO, scalapack_int N, scomplex* A, scalapack_int IA, 
          scalapack_int JA, const scalapack_desc& DESCA ) {

  scalapack_int INFO;
  pcpotri_( UPLO, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}
template <>
scalapack_int
  ppotri( const char* UPLO, scalapack_int N, dcomplex* A, scalapack_int IA, 
          scalapack_int JA, const scalapack_desc& DESCA ) {

  scalapack_int INFO;
  pzpotri_( UPLO, &N, A, &IA, &JA, DESCA.data(), &INFO );
  return INFO;

}


}
