/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemv.hpp>


SCALAPACKPP_TEST_CASE( "Gemv", "[gemv]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t M = 100, N = 200;
  const int64_t mb = 2, nb = 4;

  BlockCyclicMatrix<TestType> A( grid, M, M, mb, nb ),
                              B( grid, M, N, mb, nb ),
                              C( grid, M, N, mb, nb );

  std::fill( A.begin(), A.end(), 2 );
  std::fill( B.begin(), B.end(), 3 );
  std::fill( C.begin(), C.end(), 5 );

  pgemv( Op::NoTrans, M, M, 1.0, A.data(), 1, 1, A.desc(), 
    B.data(), 1, 1, B.desc(), 1, 2.0, C.data(), 1, 1, C.desc(), 1 );

  std::vector<TestType> C_root(M*N);
  C.gather_from( M, N, C_root.data(), M, 0, 0);
  
  const TestType val = TestType(2 * 5 + M*2*3);
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    for(int k = 0; k < N; ++k)
    for(int i = 0; i < M; ++i) {
      if( k == 0 ) CHECK( std::real(C_root[i + k*M]) == Approx(std::real(val)) );
      else         CHECK( C_root[i + k*M] == TestType(5) );
    }
  }


  MPI_Barrier(MPI_COMM_WORLD);
}
