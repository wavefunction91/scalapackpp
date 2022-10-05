/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/symm.hpp>
#include <scalapackpp/pblas/hemm.hpp>


SCALAPACKPP_REAL_TEST_CASE( "Symm", "[symm]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t M = 100, N = 200;

  BlockCyclicDist2D mat_dist( grid, 2, 4 );

  auto [M_loc1, N_loc1] = mat_dist.get_local_dims( M, M );
  auto [M_loc2, N_loc2] = mat_dist.get_local_dims( M, N );


  std::vector< TestType > A_local( M_loc1 * N_loc1, TestType(2) );
  std::vector< TestType > B_local( M_loc2 * N_loc2, TestType(3) );
  std::vector< TestType > C_local( M_loc2 * N_loc2, TestType(5) );


  auto desc_a = mat_dist.descinit_noerror( M, M, M_loc1 );
  auto desc_b = mat_dist.descinit_noerror( M, N, M_loc2 );


  psymm( 
    Side::Left, blacspp::Uplo::Lower, M, N, 
    1., A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b,
    2., C_local.data(), 1, 1, desc_b 
  );


  const TestType val = TestType(2 * 5 + M*2*3);
  for( auto x : C_local )
    CHECK( std::real(x) == Approx( std::real(val) ) );

  MPI_Barrier(MPI_COMM_WORLD);
}

SCALAPACKPP_TEST_CASE( "Hemm", "[hemm]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t M = 100, N = 200;

  BlockCyclicDist2D mat_dist( grid, 2, 4 );

  auto [M_loc1, N_loc1] = mat_dist.get_local_dims( M, M );
  auto [M_loc2, N_loc2] = mat_dist.get_local_dims( M, N );


  std::vector< TestType > A_local( M_loc1 * N_loc1, TestType(2) );
  std::vector< TestType > B_local( M_loc2 * N_loc2, TestType(3) );
  std::vector< TestType > C_local( M_loc2 * N_loc2, TestType(5) );


  auto desc_a = mat_dist.descinit_noerror( M, M, M_loc1 );
  auto desc_b = mat_dist.descinit_noerror( M, N, M_loc2 );


  phemm( 
    Side::Left, blacspp::Uplo::Lower, M, N,
    1., A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b,
    2., C_local.data(), 1, 1, desc_b 
  );


  const TestType val = TestType(2 * 5 + M*2*3);
  for( auto x : C_local )
    CHECK( std::real(x) == Approx( std::real(val) ) );

  MPI_Barrier(MPI_COMM_WORLD);
}
