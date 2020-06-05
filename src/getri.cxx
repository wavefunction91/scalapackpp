/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/matrix_inverse/getri.hpp>


using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void psgetri_( const scalapack_int* N, float * A, const scalapack_int* IA, const scalapack_int* JA, 
          const scalapack_int* DESCA, const scalapack_int* IPIV, 
          float* WORK, const scalapack_int* LWORK, scalapack_int* IWORK, const scalapack_int* LIWORK, 
          scalapack_int* INFO );
void pdgetri_( const scalapack_int* N, double * A, const scalapack_int* IA, const scalapack_int* JA, 
          const scalapack_int* DESCA, const scalapack_int* IPIV, 
          double* WORK, const scalapack_int* LWORK, scalapack_int* IWORK, const scalapack_int* LIWORK, 
          scalapack_int* INFO );
void pcgetri_( const scalapack_int* N, scomplex * A, const scalapack_int* IA, const scalapack_int* JA, 
          const scalapack_int* DESCA, const scalapack_int* IPIV, 
          scomplex* WORK, const scalapack_int* LWORK, scalapack_int* IWORK, const scalapack_int* LIWORK, 
          scalapack_int* INFO );
void pzgetri_( const scalapack_int* N, dcomplex * A, const scalapack_int* IA, const scalapack_int* JA, 
          const scalapack_int* DESCA, const scalapack_int* IPIV, 
          dcomplex* WORK, const scalapack_int* LWORK, scalapack_int* IWORK, const scalapack_int* LIWORK, 
          scalapack_int* INFO );

}


namespace scalapackpp::wrappers {

template <>
scalapack_int
  pgetri( scalapack_int N, float* A, scalapack_int IA, scalapack_int JA, 
          const scalapack_desc& DESCA, const scalapack_int* IPIV, 
          float* WORK, scalapack_int LWORK, scalapack_int* IWORK, scalapack_int LIWORK ) {

  scalapack_int INFO;
  psgetri_( &N, A, &IA, &JA, DESCA.data(), IPIV, WORK, &LWORK, IWORK, &LIWORK, &INFO );
  return INFO;

}
template <>
scalapack_int
  pgetri( scalapack_int N, double* A, scalapack_int IA, scalapack_int JA, 
	  const scalapack_desc& DESCA, const scalapack_int* IPIV, 
	  double* WORK, scalapack_int LWORK, scalapack_int* IWORK, scalapack_int LIWORK ) {

  scalapack_int INFO;
  pdgetri_( &N, A, &IA, &JA, DESCA.data(), IPIV, WORK, &LWORK, IWORK, &LIWORK, &INFO );
  return INFO;

}

template <>
scalapack_int
  pgetri( scalapack_int N, scomplex* A, scalapack_int IA, scalapack_int JA, 
          const scalapack_desc& DESCA, const scalapack_int* IPIV, 
          scomplex* WORK, scalapack_int LWORK, scalapack_int* IWORK, scalapack_int LIWORK ) {

  scalapack_int INFO;
  pcgetri_( &N, A, &IA, &JA, DESCA.data(), IPIV, WORK, &LWORK, IWORK, &LIWORK, &INFO );
  return INFO;

}

template <>
scalapack_int
  pgetri( scalapack_int N, dcomplex* A, scalapack_int IA, scalapack_int JA, 
          const scalapack_desc& DESCA, const scalapack_int* IPIV, 
          dcomplex* WORK, scalapack_int LWORK, scalapack_int* IWORK, scalapack_int LIWORK ) {

  scalapack_int INFO;
  pzgetri_( &N, A, &IA, &JA, DESCA.data(), IPIV, WORK, &LWORK, IWORK, &LIWORK, &INFO );
  return INFO;

}


}




