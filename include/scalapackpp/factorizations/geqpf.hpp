/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/factorizations/geqpf.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <scalapackpp/information.hpp>

namespace scalapackpp {

template <
  typename T,
  detail::enable_if_scalapack_real_supported_t<T, bool> = true
>
int64_t
  pgeqpf( int64_t M, int64_t N, T* A, int64_t IA, int64_t JA,
          const scalapack_desc& DESCA, int64_t* IPIV, T* TAU ) {

  auto LOCC_PIV = local_col_from_desc( JA + N - 1, DESCA );
  std::vector<internal::scalapack_int> _IPIV( LOCC_PIV );


  int64_t LWORK = -1;
  std::vector< T > WORK(5);

  wrappers::pgeqpf( M, N, A, IA, JA, DESCA, _IPIV.data(), TAU, WORK.data(), LWORK );

  LWORK = int64_t( std::real(WORK[0]) );
  WORK.resize( LWORK );

  auto INFO =  wrappers::pgeqpf( M, N, A, IA, JA, DESCA, _IPIV.data(), TAU, 
                                 WORK.data(), LWORK );

  for( int64_t i = 0; i < _IPIV.size(); ++i ) IPIV[i] = _IPIV[i];

  return INFO;

}

template <
  typename T,
  detail::enable_if_scalapack_complex_supported_t<T, bool> = true
>
int64_t
  pgeqpf( int64_t M, int64_t N, T* A, int64_t IA, int64_t JA,
          const scalapack_desc& DESCA, int64_t* IPIV, T* TAU ) {

  auto LOCC_PIV = local_col_from_desc( JA + N - 1, DESCA );
  std::vector<internal::scalapack_int> _IPIV( LOCC_PIV );


  int64_t LWORK = -1, LRWORK = -1;
  std::vector< T > WORK(5);
  std::vector< detail::real_t<T> > RWORK(5);

  wrappers::pgeqpf( M, N, A, IA, JA, DESCA, _IPIV.data(), TAU, WORK.data(), LWORK,
                    RWORK.data(), LRWORK );

  LWORK  = int64_t( std::real(WORK[0]) );
  LRWORK = int64_t( RWORK[0] );
  WORK.resize( LWORK );
  RWORK.resize( LRWORK );

  auto INFO =  wrappers::pgeqpf( M, N, A, IA, JA, DESCA, _IPIV.data(), TAU, 
                                 WORK.data(), LWORK, RWORK.data(), LRWORK );

  for( int64_t i = 0; i < _IPIV.size(); ++i ) IPIV[i] = _IPIV[i];

  return INFO;

}

}


