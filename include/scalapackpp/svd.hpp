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
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  pgesvd( Job jobu, Job jobvt, int64_t M, int64_t N,
         T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         T* S,
         T* U,  int64_t IU,  int64_t JU,  const scalapack_desc& DESCU, 
         T* VT, int64_t IVT, int64_t JVT, const scalapack_desc& DESCVT
  ) {

  auto JOBU  = char( jobu  );
  auto JOBVT = char( jobvt );

  int64_t LWORK = -1;
  std::vector< T > WORK( 5 );

  wrappers::pgesvd( &JOBU, &JOBVT, M, N, A, IA, JA, DESCA, S,
    U, IU, JU, DESCU, VT, IVT, JVT, DESCVT, WORK.data(), LWORK );

  LWORK = int64_t( WORK[0] );
  WORK.resize( LWORK );

  return wrappers::pgesvd( &JOBU, &JOBVT, M, N, A, IA, JA, DESCA, S,
    U, IU, JU, DESCU, VT, IVT, JVT, DESCVT, WORK.data(), LWORK );

}



template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, int64_t>
  pgesvd( Job jobu, Job jobvt, int64_t M, int64_t N,
         T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         detail::real_t<T>* S,
         T* U,  int64_t IU,  int64_t JU,  const scalapack_desc& DESCU, 
         T* VT, int64_t IVT, int64_t JVT, const scalapack_desc& DESCVT
  ) {

     
  auto JOBU  = char( jobu  );
  auto JOBVT = char( jobvt );

  int64_t LWORK = -1;
  std::vector< T > WORK( 5 );

  int64_t LRWORK = 1 + 4 * std::min(M,N);
  std::vector< detail::real_t<T> > RWORK( LRWORK );

  wrappers::pgesvd( &JOBU, &JOBVT, M, N, A, IA, JA, DESCA, S,
    U, IU, JU, DESCU, VT, IVT, JVT, DESCVT, WORK.data(), LWORK, RWORK.data() );

  LWORK  = int64_t(std::real(WORK[0]));
  LRWORK = int64_t(RWORK[0]); 
  WORK.resize( LWORK );
  RWORK.resize(LRWORK);

  return wrappers::pgesvd( &JOBU, &JOBVT, M, N, A, IA, JA, DESCA, S,
    U, IU, JU, DESCU, VT, IVT, JVT, DESCVT, WORK.data(), LWORK, RWORK.data() );

}


template <typename T>
detail::enable_if_scalapack_supported_t<T,int64_t>
  pgesvd( Job jobu, Job jobvt, BlockCyclicMatrix<T>& A,
          detail::real_t<T>* S, BlockCyclicMatrix<T>& U, BlockCyclicMatrix<T>& VT ) {

  return pgesvd( jobu, jobvt, A.m(), A.n(), A.data(), 1, 1, A.desc(),
                 S, U.data(), 1, 1, U.desc(), VT.data(), 1, 1, VT.desc() );

}
           
}
