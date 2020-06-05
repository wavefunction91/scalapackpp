/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/linear_systems/trtrs.hpp>
#include <scalapackpp/pblas/trsm.hpp>
#include <scalapackpp/fill_triangle.hpp>


SCALAPACKPP_TEST_CASE( "Trtrs", "[trtrs]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapack_int M = 100, N = 200;

  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc1, N_loc1] = mat_dist.get_local_dims( M, M );
  auto [M_loc2, N_loc2] = mat_dist.get_local_dims( M, N );


  std::vector< TestType > A_local( M_loc1*N_loc1, 1 );
  std::vector< TestType > B_local( M_loc2*N_loc2, 2 );
  std::vector< TestType > B_ref_local( B_local );

  // Zero out upper triangle
  fill_triangle( mat_dist, blacspp::Triangle::Upper, M, M, 
    A_local.data(), M_loc1, 0. );


  // Compute reference solution with GEMM
  auto desc_a = mat_dist.descinit_noerror( M, M, M_loc1 );
  auto desc_b = mat_dist.descinit_noerror( M, N, M_loc2 );

  // Compute reference sol with trsm
  ptrsm(
    SideFlag::Left, blacspp::Triangle::Lower, TransposeFlag::NoTranspose,
    blacspp::Diagonal::Unit, M, N, 1., A_local.data(), 1, 1, desc_a,
    B_ref_local.data(), 1, 1, desc_b
  );

  // Compute with trtrs ( B <- A**-1*B )
  auto info = ptrtrs(
    blacspp::Triangle::Lower,
    TransposeFlag::NoTranspose, blacspp::Diagonal::Unit,
    M, N, A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b 
  );

  REQUIRE( info == 0 );

  // Check B = C
  for( auto i = 0; i < B_local.size(); ++i )
    CHECK( std::real(B_local[i]) == Approx( std::real(B_ref_local[i]) ) );

}

