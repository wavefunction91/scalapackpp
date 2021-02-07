/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/linear_systems/gels.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgels( Op trans, int64_t M, int64_t N, int64_t NRHS, 
    T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
    T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {

  auto TRANS = char( trans );

  int64_t LWORK = -1;
  std::vector< T > WORK(5);

  wrappers::pgels( &TRANS, M, N, NRHS, A, IA, JA, DESCA, 
                   B, IB, JB, DESCB, WORK.data(), LWORK );

  LWORK = int64_t( std::real(WORK[0]) );
  WORK.resize( LWORK );

  return wrappers::pgels( &TRANS, M, N, NRHS, A, IA, JA, DESCA, 
                          B, IB, JB, DESCB, WORK.data(), LWORK );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  pgels( Op trans, BlockCyclicMatrix<T>& A, BlockCyclicMatrix<T>& B ) {

  // TODO sanity check
  return pgels( trans, A.m(), A.n(), B.n(), A.data(), 1, 1, A.desc(),
                B.data(), 1, 1, B.desc() );

}

}
