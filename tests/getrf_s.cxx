/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/factorizations/getrf.hpp>
#include <scalapackpp/linear_systems/getrs.hpp>







SCALAPACKPP_TEST_CASE( "Getrf_s", "[getrf_s]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100, NRHS = 20;

  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );
  auto [xxxxx, NRHS_loc ] = mat_dist.get_local_dims( M, NRHS );

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > B_local( M_loc * NRHS_loc );
  std::vector< TestType > A_SPD_local( M_loc * N_loc );

  std::generate( A_local.begin(), A_local.end(), [&](){ return dist(gen); } );
  std::generate( B_local.begin(), B_local.end(), [&](){ return dist(gen); } );
  auto B_copy( B_local );

  for( auto i = 0; i < M; ++i ) 
  if( mat_dist.i_own( i, i ) ) {

    auto [ row_idx, col_idx ] = mat_dist.local_indx(i, i);
    A_local[ row_idx + col_idx * M_loc ] += 10.;

  }

  auto A_copy( A_local );

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );
  auto desc_rhs = mat_dist.descinit_noerror( M, NRHS, M_loc );


  // Solve linear system
  std::vector<int64_t> IPIV( M_loc + mat_dist.mb() );
  auto info = pgetrf( M, M, A_local.data(), 1, 1, desc, IPIV.data() );
  REQUIRE( info == 0 );

  info = pgetrs( TransposeFlag::NoTranspose, M, NRHS, A_local.data(), 1, 1, desc, IPIV.data(),
            B_local.data(), 1, 1, desc_rhs );

  REQUIRE( info == 0 );
  // Check correctness
  pgemm( TransposeFlag::NoTranspose, TransposeFlag::NoTranspose, M, NRHS, M,
         -1., A_copy.data(), 1, 1, desc, B_local.data(), 1, 1, desc_rhs,
         1. , B_copy.data(), 1, 1, desc_rhs );

  auto tol = 100*M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
  for( auto x : B_copy ) CHECK( std::abs(x) < tol );

}
