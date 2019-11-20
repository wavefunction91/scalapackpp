#include <scalapackpp/wrappers/factorizations/getrf.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void psgetrf_( const scalapack_int* M, const scalapack_int* N, float* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, scalapack_int* INFO );
void pdgetrf_( const scalapack_int* M, const scalapack_int* N, double* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, scalapack_int* INFO );
void pcgetrf_( const scalapack_int* M, const scalapack_int* N, scomplex* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, scalapack_int* INFO );
void pzgetrf_( const scalapack_int* M, const scalapack_int* N, dcomplex* A, 
  const scalapack_int* IA, const scalapack_int* JA, const scalapack_int* DESCA, 
  scalapack_int* IPIV, scalapack_int* INFO );

}


namespace scalapackpp::wrappers {

template <>
scalapack_int
  pgetrf( scalapack_int M, scalapack_int N, float* A, scalapack_int IA, scalapack_int JA,
          const scalapack_desc& DESCA, scalapack_int* IPIV ) {

  scalapack_int INFO;
  psgetrf_( &M, &N, A, &IA, &JA, DESCA.data(), IPIV, &INFO );
  return INFO;

}

template <>
scalapack_int
  pgetrf( scalapack_int M, scalapack_int N, double* A, scalapack_int IA, scalapack_int JA,
          const scalapack_desc& DESCA, scalapack_int* IPIV ) {

  scalapack_int INFO;
  pdgetrf_( &M, &N, A, &IA, &JA, DESCA.data(), IPIV, &INFO );
  return INFO;

}

template <>
scalapack_int
  pgetrf( scalapack_int M, scalapack_int N, scomplex* A, scalapack_int IA, scalapack_int JA,
          const scalapack_desc& DESCA, scalapack_int* IPIV ) {

  scalapack_int INFO;
  pcgetrf_( &M, &N, A, &IA, &JA, DESCA.data(), IPIV, &INFO );
  return INFO;

}

template <>
scalapack_int
  pgetrf( scalapack_int M, scalapack_int N, dcomplex* A, scalapack_int IA, scalapack_int JA,
          const scalapack_desc& DESCA, scalapack_int* IPIV ) {

  scalapack_int INFO;
  pzgetrf_( &M, &N, A, &IA, &JA, DESCA.data(), IPIV, &INFO );
  return INFO;

}

}
