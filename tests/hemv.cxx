/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/symv.hpp>
#include <scalapackpp/pblas/hemv.hpp>


SCALAPACKPP_REAL_TEST_CASE( "Symv", "[symv]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t N = 200;
  const int64_t mb = 2, nb = 4;

  BlockCyclicMatrix<TestType> A( grid, N, N, mb, nb ),
                              B( grid, N, N, mb, nb ),
                              C( grid, N, N, mb, nb );

  std::fill( A.begin(), A.end(), 2 );
  std::fill( B.begin(), B.end(), 3 );
  std::fill( C.begin(), C.end(), 5 );

  psymv( Uplo::Lower, N, 1.0, A.data(), 1, 1, A.desc(), 
    B.data(), 1, 1, B.desc(), 1, 2.0, C.data(), 1, 1, C.desc(), 1 );

  std::vector<TestType> C_root(N*N);
  C.gather_from( N, N, C_root.data(), N, 0, 0);
  
  const TestType val = TestType(2 * 5 + N*2*3);
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    for(int k = 0; k < N; ++k)
    for(int i = 0; i < N; ++i) {
      if( k == 0 ) CHECK( std::real(C_root[i + k*N]) == Approx(std::real(val)) );
      else         CHECK( C_root[i + k*N] == TestType(5) );
    }
  }


  MPI_Barrier(MPI_COMM_WORLD);
}

SCALAPACKPP_REAL_TEST_CASE( "Hemv", "[hemv]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t N = 200;
  const int64_t mb = 2, nb = 4;

  BlockCyclicMatrix<TestType> A( grid, N, N, mb, nb ),
                              B( grid, N, N, mb, nb ),
                              C( grid, N, N, mb, nb );

  std::fill( A.begin(), A.end(), 2 );
  std::fill( B.begin(), B.end(), 3 );
  std::fill( C.begin(), C.end(), 5 );

  phemv( Uplo::Lower, N, 1.0, A.data(), 1, 1, A.desc(), 
    B.data(), 1, 1, B.desc(), 1, 2.0, C.data(), 1, 1, C.desc(), 1 );

  std::vector<TestType> C_root(N*N);
  C.gather_from( N, N, C_root.data(), N, 0, 0);
  
  const TestType val = TestType(2 * 5 + N*2*3);
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    for(int k = 0; k < N; ++k)
    for(int i = 0; i < N; ++i) {
      if( k == 0 ) CHECK( std::real(C_root[i + k*N]) == Approx(std::real(val)) );
      else         CHECK( C_root[i + k*N] == TestType(5) );
    }
  }


  MPI_Barrier(MPI_COMM_WORLD);
}
