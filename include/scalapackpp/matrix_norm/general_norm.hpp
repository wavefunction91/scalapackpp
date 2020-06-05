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

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, detail::real_t<T>>
  general_norm( const BlockCyclicDist2D& dist, MatrixNorm norm, 
                scalapack_int M, scalapack_int N, const T* A, scalapack_int IA, 
                scalapack_int JA, const scalapack_desc& DESCA ) {

  auto [ ip, jp ]   = dist.owner_coordinate( IA-1, JA-1 );

  auto iro = (IA-1) % dist.mb();
  auto ico = (JA-1) % dist.nb();

  auto Mp0 = numroc( M+iro, dist.mb(), dist.grid().ipr(), ip, dist.grid().npr() );
  auto Nq0 = numroc( N+ico, dist.nb(), dist.grid().ipc(), jp, dist.grid().npc() );

  scalapack_int LWORK = 0;
  if( norm == MatrixNorm::OneNorm )           LWORK = Nq0;
  else if( norm == MatrixNorm::InfinityNorm ) LWORK = Mp0;

  LWORK = std::max( 1, LWORK );
  std::vector< detail::real_t<T> > WORK( LWORK );

  auto NORM = detail::type_string( norm );
  return wrappers::plange( NORM.c_str(), M, N, A, IA, JA, DESCA, WORK.data() );

}

}

