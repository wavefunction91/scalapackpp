/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/matrix_inverse/getri.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  pgetri( scalapack_int N, T* A, scalapack_int IA, scalapack_int JA, 
          const scalapack_desc& DESCA, scalapack_int* IPIV ) { 


  scalapack_int LWORK = -1, LIWORK = -1;
  std::vector< T > WORK( 5 );
  std::vector< scalapack_int > IWORK(5);

  auto info = wrappers::pgetri( N, A, IA, JA, DESCA, IPIV, WORK.data(), 
    LWORK, IWORK.data(), LIWORK );

  LWORK = scalapack_int( std::real(WORK[0]) );
  LIWORK = IWORK[0];

  if( LWORK > 0 and LIWORK > 0 ) {
    WORK.resize(LWORK);
    IWORK.resize(LIWORK);
    info = wrappers::pgetri( N, A, IA, JA, DESCA, IPIV, WORK.data(), 
      LWORK, IWORK.data(), LIWORK );
  }

  return info;
}


}




