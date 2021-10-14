/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/lacpy.hpp>

SCALAPACKPP_TEST_CASE( "Lacpy", "[lacpy]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t M = 100, N = 200;
  const int64_t mb = 2, nb = 4;

  BlockCyclicMatrix<TestType> A_sca( grid, M, N, mb, nb ),
                              B_sca( grid, M, N, mb, nb );

  std::vector< TestType > A;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    A.resize( M*N, 0 );

    for( int j = 0; j < N; ++j )
    for( int i = 0; i < M; ++i )
      if( i == j )      A[ i + j*M ] = 1.;
      else if( i < j )  A[ i + j*M ] = 2.;
      else              A[ i + j*M ] = 3.;

  }

  A_sca.scatter_to( M, N, A.data(), M, 0, 0 );

  std::fill( B_sca.begin(), B_sca.end(), 0. );
  placpy( Uplo::Upper, A_sca, B_sca );

  // Gather results to A
  B_sca.gather_from( M, N, A.data(), M, 0, 0 );

  // Check
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    for( int j = 0; j < N; ++j )
    for( int i = 0; i < M; ++i ) {
      if( i == j )      CHECK( std::real(A[ i + j*M ]) == Approx(1.) );
      else if( i < j )  CHECK( std::real(A[ i + j*M ]) == Approx(2.) );
      else              CHECK( std::real(A[ i + j*M ]) == Approx(0.) );
    }

  }
}

