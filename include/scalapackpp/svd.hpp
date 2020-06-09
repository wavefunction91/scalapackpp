/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/gesvd.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, scalapack_int>
  pgesvd( VectorFlag jobu, VectorFlag jobvt, scalapack_int M, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         T* S,
         T* U,  scalapack_int IU,  scalapack_int JU,  const scalapack_desc& DESCU, 
         T* VT, scalapack_int IVT, scalapack_int JVT, const scalapack_desc& DESCVT
  ) {

  auto JOBU  = detail::type_string( jobu  );
  auto JOBVT = detail::type_string( jobvt );

  scalapack_int LWORK = -1;
  std::vector< T > WORK( 5 );

  wrappers::pgesvd( JOBU.c_str(), JOBVT.c_str(), M, N, A, IA, JA, DESCA, S,
    U, IU, JU, DESCU, VT, IVT, JVT, DESCVT, WORK.data(), LWORK );

  LWORK = scalapack_int( WORK[0] );
  WORK.resize( LWORK );

  return wrappers::pgesvd( JOBU.c_str(), JOBVT.c_str(), M, N, A, IA, JA, DESCA, S,
    U, IU, JU, DESCU, VT, IVT, JVT, DESCVT, WORK.data(), LWORK );

}



template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, scalapack_int>
  pgesvd( VectorFlag jobu, VectorFlag jobvt, scalapack_int M, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         detail::real_t<T>* S,
         T* U,  scalapack_int IU,  scalapack_int JU,  const scalapack_desc& DESCU, 
         T* VT, scalapack_int IVT, scalapack_int JVT, const scalapack_desc& DESCVT
  ) {

     
  auto JOBU  = detail::type_string( jobu  );
  auto JOBVT = detail::type_string( jobvt );

  scalapack_int LWORK = -1;
  std::vector< T > WORK( 5 );

  scalapack_int LRWORK = 1 + 4 * std::min(M,N);
  std::vector< detail::real_t<T> > RWORK( LRWORK );

  wrappers::pgesvd( JOBU.c_str(), JOBVT.c_str(), M, N, A, IA, JA, DESCA, S,
    U, IU, JU, DESCU, VT, IVT, JVT, DESCVT, WORK.data(), LWORK, RWORK.data() );

  LWORK  = scalapack_int(std::real(WORK[0]));
  LRWORK = scalapack_int(RWORK[0]); 
  WORK.resize( LWORK );
  RWORK.resize(LRWORK);

  return wrappers::pgesvd( JOBU.c_str(), JOBVT.c_str(), M, N, A, IA, JA, DESCA, S,
    U, IU, JU, DESCU, VT, IVT, JVT, DESCVT, WORK.data(), LWORK, RWORK.data() );

}
           
}
