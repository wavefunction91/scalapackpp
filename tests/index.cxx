/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/scatter_gather.hpp>

TEST_CASE( "Index Conversion", "[index]" ) {


  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100;

  BlockCyclicDist2D mat_dist( grid, 4, 2 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );

  SECTION( "Local Index" ) {

    std::vector<double> A_local( M_loc*N_loc, 0. );

    // Create distributed identity
    for( auto i = 0; i < M; ++i ) 
    if( mat_dist.i_own( i, i ) ) {

      auto [ row_idx, col_idx ] = mat_dist.local_indx(i, i);
      A_local[ row_idx + col_idx * M_loc ] = 1.;

    }

    std::vector<double> A;
    if( grid.ipr() == 0 and grid.ipc() == 0 ) 
      A.resize( M*M, 0. );

    mat_dist.gather( M, M, A.data(), M, A_local.data(), M_loc, 0, 0 );
    
    if( grid.ipr() == 0 and grid.ipc() == 0 ) {
      for( auto i = 0; i < M; ++i )
      for( auto j = 0; j < M; ++j ) {
        if( i == j ) CHECK( A[i + j*M] == 1. );
        else         CHECK( A[i + j*M] == 0. );
      }
    }

  }


  MPI_Barrier(MPI_COMM_WORLD);
}
