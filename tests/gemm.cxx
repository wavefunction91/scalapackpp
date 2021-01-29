/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>


SCALAPACKPP_TEST_CASE( "Gemm", "[gemm]" ) {

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

  pgemm( TransposeFlag::NoTranspose, TransposeFlag::NoTranspose,
         1., A, B, 2., C );

  const TestType val = TestType(2 * 5 + M*2*3);
  for( auto x : C )
    CHECK( std::real(x) == Approx( std::real(val) ) );


}
