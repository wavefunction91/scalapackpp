/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/information.hpp>

namespace scalapackpp {

int64_t numroc( int64_t N, int64_t NB,
  int64_t IPROC, int64_t ISRC, 
  int64_t NPROC ) {

  auto dist    = ( NPROC + IPROC - ISRC ) % NPROC;
  auto nblocks = N / NB;
  auto dim     = ( nblocks / NPROC ) * NB;
  auto extra   = nblocks % NPROC;

  if( dist < extra )       dim += NB;     // Extra block
  else if( dist == extra ) dim += N % NB; // Last block

  return dim;
}

std::pair<int64_t, int64_t>
  get_local_dims( const blacspp::Grid& grid,
    int64_t M, int64_t N,
    int64_t MB, int64_t NB,
    int64_t ISRC, int64_t JSRC ) {

  return { numroc( M, MB, grid.ipr(), ISRC, grid.npr() ),
           numroc( N, NB, grid.ipc(), JSRC, grid.npc() ) };

}

/*
std::pair< block_cyclic_coordinate, block_cyclic_coordinate >
  local_index( const blacspp::Grid& grid,
    scalapack_int MB, scalapack_int NB,
    scalapack_int I, scalapack_int J
  ) {

  block_cyclic_coordinate row_coord, col_coord;

  row_coord.block_idx = I / ( grid.npr() * MB );
  col_coord.block_idx = J / ( grid.npc() * NB );

  row_coord.process_id = (I / MB) & grid.npr();
  col_coord.process_id = (J / NB) & grid.npc();

  row_coord.local_idx  = I % MB;
  col_coord.local_idx  = J % NB;

  return {row_coord, col_coord};

}
*/



}
