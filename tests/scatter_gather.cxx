/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <blacspp/information.hpp>


SCALAPACKPP_TEST_CASE( "Gather", "[scatter-gather]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t M = 100, N = 200;
  BlockCyclicDist2D mat_dist( grid, 2, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, N );

  std::vector< TestType > data_local( M_loc * N_loc, TestType(mpi.rank()) );
  std::vector< TestType > data_gathered;
  if( grid.ipr() == 0 and grid.ipc() == 0 )
    data_gathered.resize( M*N, TestType(-1) );


  mat_dist.gather( M, N, data_gathered.data(), M, data_local.data(), M_loc, 0, 0 );

  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    for( int i = 0; i < M; ++i )
    for( int j = 0; j < N; ++j ) {

      auto rank = blacspp::coordinate_rank( grid, mat_dist.owner_coordinate(i,j) );
      CHECK( data_gathered[i + j*M] == TestType(rank) );

    }

  } else {
    CHECK( data_gathered.size() == 0 );
  }

  for( auto x : data_local ) CHECK( x == TestType(mpi.rank()) );

  grid.barrier( blacspp::Scope::All );

  MPI_Barrier(MPI_COMM_WORLD);
}



SCALAPACKPP_TEST_CASE( "Scatter", "[scatter-gather]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t M = 100, N = 200;
  BlockCyclicDist2D mat_dist( grid, 2, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, N );


  std::vector< TestType > data_local( M_loc * N_loc, TestType(-1) );
  std::vector< TestType > data_gathered;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    data_gathered.resize( M*N );

    for( int i = 0; i < M; ++i )
    for( int j = 0; j < N; ++j ) {
      data_gathered[i + j*M] = 
        blacspp::coordinate_rank( grid, mat_dist.owner_coordinate(i,j) );
    }

  }


  mat_dist.scatter( M, N, data_gathered.data(), M, data_local.data(), M_loc, 0, 0 );

  for( auto x : data_local ) CHECK( x == TestType( mpi.rank() ) );

  MPI_Barrier(MPI_COMM_WORLD);
}
