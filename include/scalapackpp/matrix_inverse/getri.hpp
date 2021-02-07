/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/matrix_inverse/getri.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgetri( int64_t N, T* A, int64_t IA, int64_t JA, 
          const scalapack_desc& DESCA, int64_t* IPIV ) { 

  auto LOCR_A = local_row_from_desc( DESCA[internal::_M_A], DESCA );

  std::vector<internal::scalapack_int> _IPIV( LOCR_A + DESCA[internal::_MB_A] );
  for( int64_t i = 0; i < _IPIV.size(); ++i )
    _IPIV[i] = detail::to_scalapack_int(IPIV[i]) ;

  int64_t LWORK = -1, LIWORK = -1;
  std::vector< T > WORK( 5 );
  std::vector< internal::scalapack_int > IWORK(5);

  auto info = wrappers::pgetri( N, A, IA, JA, DESCA, _IPIV.data(), WORK.data(), 
    LWORK, IWORK.data(), LIWORK );

  LWORK = int64_t( std::real(WORK[0]) );
  LIWORK = IWORK[0];

  if( LWORK > 0 and LIWORK > 0 ) {
    WORK.resize(LWORK);
    IWORK.resize(LIWORK);
    info = wrappers::pgetri( N, A, IA, JA, DESCA, _IPIV.data(), WORK.data(), 
      LWORK, IWORK.data(), LIWORK );
  }

  return info;
}



template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgetri( BlockCyclicMatrix<T>& A, int64_t* IPIV ) { 

  // TODO sanity check
  return pgetri( A.m(), A.data(), 1, 1, A.desc(), IPIV );

}


}




