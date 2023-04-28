/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>
#include <scalapackpp/pblas/copy.hpp>

SCALAPACKPP_TEST_CASE( "Copy", "[copy]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t M = 100, N = 200;
  const int64_t mb = 2, nb = 4;

  BlockCyclicMatrix<TestType> A_sca( grid, M, N, mb, nb );
  BlockCyclicMatrix<TestType> C_sca( grid, M, N, mb, nb );

  std::vector< TestType > A, B;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    A.resize( M*N, 0 );
    B.resize( M*N, 0 );

    for( int j = 0; j < N; ++j )
    for( int i = 0; i < M; ++i ) {
      if( i == j )      A[ i + j*M ] = 1.;
      else if( i < j )  A[ i + j*M ] = 2.;
      else              A[ i + j*M ] = 3.;
    }

  }

  for(int j = 0; j < N; ++j) {
    A_sca.scatter_to( M, N, A.data(), M, 0, 0 );
    C_sca.scatter_to( M, N, A.data(), M, 0, 0 );

    pcopy(M, C_sca.data(), 1, 1, C_sca.desc(), 1, A_sca.data(), 1, j+1, A_sca.desc(), 1);

    // Gather results to B
    A_sca.gather_from( M, N, B.data(), M, 0, 0 );

    // Check
    if( grid.ipr() == 0 and grid.ipc() == 0 ) {

      for(int k = 0; k < N; ++k)
      for(int i = 0; i < M; ++i) {
        if(k == j) CHECK( B[i + k*M] == A[i] );
        else       CHECK( A[i + k*M] == B[i + k*M] );
      }

    }
        
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

