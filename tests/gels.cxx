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
  int64_t mb = 4;

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A,B;

  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    A.resize( M*N ); B.resize( M*NRHS );
    std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );
    std::generate( B.begin(), B.end(), [&](){ return dist(gen); } );

  }

  BlockCyclicMatrix<TestType> A_sca( grid, M, N,    mb, mb ),
                              B_sca( grid, M, NRHS, mb, mb );

  A_sca.scatter_to( M, N,    A.data(), M, 0, 0 );
  B_sca.scatter_to( M, NRHS, B.data(), M, 0, 0 );

  // Solve linear system
  auto info = pgels( Op::NoTrans, A_sca, B_sca );
  REQUIRE( info == 0 );

  // Check correctness
  std::vector< TestType > X;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    X.resize( M * NRHS );
  }

  B_sca.gather_from( M, NRHS, X.data(), M, 0, 0 );

  auto tol = M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    lapack::gels( lapack::Op::NoTrans, M, N, NRHS, A.data(), M, B.data(), M );
    for( int i = 0; i < M*NRHS; ++i ) 
      CHECK( std::abs(X[i] - B[i]) < tol );
  }

  MPI_Barrier(MPI_COMM_WORLD);
}
