/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/factorizations/geqrf.hpp>
#include <scalapackpp/factorizations/geqpf.hpp>
#include <scalapackpp/householder/generate_q_householder.hpp>
#include <scalapackpp/pblas/gemm.hpp>

SCALAPACKPP_TEST_CASE( "Geqrf", "[geqrf]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100, N = 50;


  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, N );
  auto [RM_loc, RN_loc] = mat_dist.get_local_dims( N, N );
  auto desc_a = mat_dist.descinit_noerror( M, N, M_loc );
  auto desc_r = mat_dist.descinit_noerror( N, N, RM_loc );

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A;

  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A.resize( M*N );
    std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );
  }

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > TAU_local( local_col_from_desc( std::min(M,N), desc_a ) );

  mat_dist.scatter( M, N, A.data(), M, A_local.data(), M_loc, 0, 0 );
  std::vector< TestType > A_cpy_local(A_local);

  auto info = pgeqrf( M, N, A_local.data(), 1, 1, desc_a, TAU_local.data() );
  REQUIRE( info == 0 );

  // Generate Q
  std::vector<TestType> Q_local(A_local);
  info = generate_q_householder( M, N, N, Q_local.data(), 1, 1, desc_a, TAU_local.data() );
  REQUIRE( info == 0 );

  // Compute R = Q**H * A
  std::vector<TestType> R_local( RM_loc * RN_loc );
  pgemm( TransposeFlag::ConjTranspose, TransposeFlag::NoTranspose, N, N, M,
         TestType(1.), Q_local.data(), 1, 1, desc_a, A_cpy_local.data(), 1, 1, 
         desc_a, TestType(0.), R_local.data(), 1, 1, desc_r );



  // Check correctnesss
  std::vector<TestType> R;
  if( grid.ipc() == 0 and grid.ipr() == 0 ) {
    R.resize( N*N );
  }

  mat_dist.gather( N, N, R.data(), N, R_local.data(), RM_loc, 0, 0 );
  mat_dist.gather( M, N, A.data(), M, A_local.data(), M_loc,  0, 0 );

  if( grid.ipc() == 0 and grid.ipr() == 0 ) {

    const auto tol = M*N*std::numeric_limits<detail::real_t<TestType>>::epsilon();
    for( int j = 0; j <  N; ++j )
    for( int i = 0; i <= j; ++i ) {
      CHECK( std::abs(R[i + j*N] - A[ i + j*M ]) < tol );
    }

  }


}


SCALAPACKPP_TEST_CASE( "Geqpf", "[geqpf]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100, N = 50;


  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, N );
  auto desc_a = mat_dist.descinit_noerror( M, N, M_loc );

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A;

  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A.resize( M*N );
    std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );
  }

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > TAU_local( local_col_from_desc( std::min(M,N), desc_a ) );
  std::vector< int64_t  > IPIV_local( local_col_from_desc( N, desc_a ) );

  mat_dist.scatter( M, N, A.data(), M, A_local.data(), M_loc, 0, 0 );

  auto info = pgeqpf( M, N, A_local.data(), 1, 1, desc_a, IPIV_local.data(), 
                      TAU_local.data() );

  REQUIRE( info == 0 );

  // TODO: check for correctness
}
