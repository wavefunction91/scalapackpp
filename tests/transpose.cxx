/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/transpose.hpp>

SCALAPACKPP_TEST_CASE( "Transpose", "[trans]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t M = 100, N = 200;
  const int64_t mb = 2, nb = 4;

  BlockCyclicMatrix<TestType> A_sca( grid, M, N, mb, nb ),
                              B_sca( grid, N, M, mb, nb );


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

  SECTION( "Transpose" ) {

    transpose( A_sca, B_sca );

    // Gather results to A
    B_sca.gather_from( N, M, A.data(), N, 0, 0 );

    // Check
    if( grid.ipr() == 0 and grid.ipc() == 0 ) {


      for( int j = 0; j < M; ++j )
      for( int i = 0; i < N; ++i )
        if( i == j )      CHECK( std::real(A[ i + j*N ]) == Approx(1.) );
        else if( i > j )  CHECK( std::real(A[ i + j*N ]) == Approx(2.) );
        else              CHECK( std::real(A[ i + j*N ]) == Approx(3.) );

    }
  }

  SECTION( "ConjTranspose" ) {

    conj_transpose( A_sca, B_sca );

    // Gather results to A
    B_sca.gather_from( N, M, A.data(), N, 0, 0 );

    // Check
    if( grid.ipr() == 0 and grid.ipc() == 0 ) {


      for( int j = 0; j < M; ++j )
      for( int i = 0; i < N; ++i )
        if( i == j )      CHECK( std::real(A[ i + j*N ]) == Approx(1.) );
        else if( i > j )  CHECK( std::real(A[ i + j*N ]) == Approx(2.) );
        else              CHECK( std::real(A[ i + j*N ]) == Approx(3.) );

    }

  }
        
}

