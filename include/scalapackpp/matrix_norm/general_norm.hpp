/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/matrix_norm/lange.hpp>
#include <scalapackpp/types.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, detail::real_t<T>>
  general_norm( const BlockCyclicDist2D& dist, Norm norm, 
                int64_t M, int64_t N, const T* A, int64_t IA, 
                int64_t JA, const scalapack_desc& DESCA ) {

  auto [ ip, jp ]   = dist.owner_coordinate( IA-1, JA-1 );

  auto iro = (IA-1) % dist.mb();
  auto ico = (JA-1) % dist.nb();

  auto Mp0 = numroc( M+iro, dist.mb(), dist.grid().ipr(), ip, dist.grid().npr() );
  auto Nq0 = numroc( N+ico, dist.nb(), dist.grid().ipc(), jp, dist.grid().npc() );

  int64_t LWORK = 0;
  if( norm == Norm::One )      LWORK = Nq0;
  else if( norm == Norm::Inf ) LWORK = Mp0;

  LWORK = std::max( (int64_t)1, LWORK );
  std::vector< detail::real_t<T> > WORK( LWORK );

  auto NORM = char( norm );
  return wrappers::plange( &NORM, M, N, A, IA, JA, DESCA, WORK.data() );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T, detail::real_t<T>>
  general_norm( Norm norm, const BlockCyclicMatrix<T>& A ) {

  return general_norm( A.dist(), norm, A.m(), A.n(), A.data(), 1, 1, A.desc() );  

}

}

