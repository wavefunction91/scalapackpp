/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/information.hpp>
#include <blacspp/wrappers/support.hpp>

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


int64_t local_row_from_desc( int64_t M, const scalapack_desc& desc ) {
  const auto ICXT_A = desc[internal::_CTXT_A];
  const auto M_A    = desc[internal::_M_A];
  const auto MB_A   = desc[internal::_MB_A];
  const auto RSRC_A = desc[internal::_RSRC_A];

  auto grid_dim = blacspp::wrappers::grid_info( ICXT_A );
  return numroc( M, MB_A, grid_dim.my_row, RSRC_A, grid_dim.np_row );
}

int64_t local_col_from_desc( int64_t N, const scalapack_desc& desc ) {
  const auto ICXT_A = desc[internal::_CTXT_A];
  const auto N_A    = desc[internal::_N_A];
  const auto NB_A   = desc[internal::_NB_A];
  const auto CSRC_A = desc[internal::_CSRC_A];

  auto grid_dim = blacspp::wrappers::grid_info( ICXT_A );
  return numroc( N, NB_A, grid_dim.my_col, CSRC_A, grid_dim.np_col );
}


}
