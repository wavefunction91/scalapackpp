/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/matrix_norm/lansy.hpp>
#include <scalapackpp/types.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>
#include <scalapackpp/util/math.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, detail::real_t<T>>
  symmetric_norm( const BlockCyclicDist2D& dist, 
                  MatrixNorm norm, blacspp::Triangle uplo,
                  int64_t N, const T* A, int64_t IA, 
                  int64_t JA, const scalapack_desc& DESCA ) {

  auto [ ip, jp ]   = dist.owner_coordinate( IA-1, JA-1 );

  auto iro = (IA-1) % dist.mb();
  auto ico = (JA-1) % dist.nb();

  auto Np0 = numroc( N+iro, dist.mb(), dist.grid().ipr(), ip, dist.grid().npr() );
  auto Nq0 = numroc( N+ico, dist.nb(), dist.grid().ipc(), jp, dist.grid().npc() );

  int64_t LDW = 0;
  if( dist.grid().npr() != dist.grid().npc() ) {
    int64_t lcm = std::lcm( dist.grid().npr(), dist.grid().npc() );
    int64_t ldw = dist.mb() * detail::div_ceil(
      detail::div_ceil( Np0, dist.mb() ), (lcm/dist.grid().npr())
    );
  }

  int64_t LWORK = 0;
  if( norm == OneNorm or norm == InfinityNorm )
    LWORK = 2*(Nq0 + Np0 + LDW);

  LWORK = std::max( (int64_t)1, LWORK );
  std::vector< detail::real_t<T> > WORK( LWORK );

  auto NORM = detail::type_string( norm );
  auto UPLO = blacspp::detail::type_string( uplo );
  return wrappers::plansy( NORM.c_str(), UPLO.c_str(), N, A, IA, JA, 
    DESCA, WORK.data() );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T, detail::real_t<T>>
  symmetric_norm( MatrixNorm norm, blacspp::Triangle uplo, 
                  const BlockCyclicMatrix<T>& A ) {

  return symmetric_norm( A.dist(), norm, uplo, A.m(), A.n(), A.data(), 1, 1, 
                         A.desc() );  

}

}

