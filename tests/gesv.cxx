/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/linear_systems/gesv.hpp>

SCALAPACKPP_TEST_CASE( "Gesv", "[gesv]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100, NRHS = 20;
  int64_t mb = 4;

  std::default_random_engine gen( mpi.rank() );
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  BlockCyclicMatrix<TestType> A( grid, M, M,    mb, mb ),
                              B( grid, M, NRHS, mb, mb );

  std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );
  std::generate( B.begin(), B.end(), [&](){ return dist(gen); } );
  auto B_copy( B );

  for( auto i = 0; i < M; ++i ) 
  if( A.dist().i_own( i, i ) ) {

    auto [ row_idx, col_idx ] = A.dist().local_indx(i, i);
    A.data()[ row_idx + col_idx * A.m_local() ] += 10.;

  }

  auto A_copy( A );

  // Solve linear system
  std::vector<int64_t> IPIV( A.m_local() + A.dist().mb() );
  auto info = pgesv( A, IPIV.data(), B );

  REQUIRE( info == 0 );
  // Check correctness
  pgemm( Op::NoTrans, Op::NoTrans, 
         -1., A_copy, B, 1., B_copy );

  auto tol = 100*M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
  for( auto x : B_copy ) CHECK( std::abs(x) < tol );

  MPI_Barrier(MPI_COMM_WORLD);
}
