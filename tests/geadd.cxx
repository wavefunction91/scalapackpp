/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/geadd.hpp>

SCALAPACKPP_TEST_CASE( "Geadd", "[geadd]" ) {


  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapack_int M = 100, N = 200;

  BlockCyclicDist2D mat_dist( grid, 2, 4 );


  auto [M_loc1, N_loc1] = mat_dist.get_local_dims( M, N ); // Untranspose
  auto [M_loc2, N_loc2] = mat_dist.get_local_dims( N, M ); // Transpose


  auto desc_a = mat_dist.descinit_noerror( M, N, M_loc1 );
  auto desc_b = mat_dist.descinit_noerror( N, M, M_loc2 );


  std::vector< TestType > A_local( M_loc1 * N_loc1 );
  std::vector< TestType > B_local( M_loc2 * N_loc2 );


  std::vector< TestType > A;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    A.resize( M*N, 0 );

    for( int j = 0; j < N; ++j )
    for( int i = 0; i < M; ++i )
      if( i == j )      A[ i + j*M ] = 1.;
      else if( i < j )  A[ i + j*M ] = 2.;
      else              A[ i + j*M ] = 3.;

  }

  mat_dist.scatter( M, N, A.data(), M, A_local.data(), M_loc1, 0, 0 );

  detail::real_t<TestType> fact = 0.5;
  pgeadd( TransposeFlag::Transpose, N, M, fact, A_local.data(), 1, 1, desc_a,
          1., B_local.data(), 1, 1, desc_b );


  // Gather results to A
  mat_dist.gather( N, M, A.data(), N, B_local.data(), M_loc2, 0, 0 );

  // Check
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {


    for( int j = 0; j < M; ++j )
    for( int i = 0; i < N; ++i )
      if( i == j )      CHECK( std::real(A[ i + j*N ]) == Approx(1.*fact) );
      else if( i > j )  CHECK( std::real(A[ i + j*N ]) == Approx(2.*fact) );
      else              CHECK( std::real(A[ i + j*N ]) == Approx(3.*fact) );

  }
  
}
