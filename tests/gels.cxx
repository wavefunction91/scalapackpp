/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/linear_systems/gels.hpp>


#include <lapack.hh>

SCALAPACKPP_TEST_CASE( "Gels", "[gels]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100, N = 50, NRHS = 20;


  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc]     = mat_dist.get_local_dims( M, N );
  auto [xxxxx, NRHS_loc ] = mat_dist.get_local_dims( M, NRHS );

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A,B;

  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    A.resize( M*N ); B.resize( M*NRHS );
    std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );
    std::generate( B.begin(), B.end(), [&](){ return dist(gen); } );

  }

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > B_local( M_loc * NRHS_loc );

  mat_dist.scatter( M, N,    A.data(), M, A_local.data(), M_loc, 0, 0 );
  mat_dist.scatter( M, NRHS, B.data(), M, B_local.data(), M_loc, 0, 0 );

  // Solve linear system
  auto desc_a = mat_dist.descinit_noerror( M, N,    M_loc );
  auto desc_b = mat_dist.descinit_noerror( M, NRHS, M_loc );
  auto info = pgels( TransposeFlag::NoTranspose, M, N, NRHS, 
                     A_local.data(), 1, 1, desc_a,
                     B_local.data(), 1, 1, desc_b );

  REQUIRE( info == 0 );
  // Check correctness
  std::vector< TestType > X;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    X.resize( M * NRHS );
  }

  mat_dist.gather( M, NRHS, X.data(), M, B_local.data(), M_loc, 0, 0 );

  auto tol = M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    lapack::gels( lapack::Op::NoTrans, M, N, NRHS, A.data(), M, B.data(), M );
    for( int i = 0; i < M*NRHS; ++i ) 
      CHECK( std::abs(X[i] - B[i]) < tol );
  }

}
