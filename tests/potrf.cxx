/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/trsm.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/factorizations/potrf.hpp>







SCALAPACKPP_TEST_CASE( "Potrf", "[potrf]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int M = 100;

  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > A_SPD_local( M_loc * N_loc );

  std::generate( A_local.begin(), A_local.end(), [&](){ return dist(gen); } );

  for( auto i = 0; i < M; ++i ) 
  if( mat_dist.i_own( i, i ) ) {

    auto [ row_idx, col_idx ] = mat_dist.local_indx(i, i);
    A_local[ row_idx + col_idx * M_loc ] += 10.;

  }

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );

  // Make A SPD
  pgemm(
    TransposeFlag::ConjTranspose, TransposeFlag::NoTranspose, M, M, M, 
    1., A_local.data(), 1, 1, desc, A_local.data(), 1, 1, desc,
    0., A_SPD_local.data(), 1, 1, desc
  );
  std::vector< TestType > A_SPD_copy( A_SPD_local );


  // Perform POTRF
  auto info = ppotrf( blacspp::Triangle::Lower, M, A_SPD_local.data(), 1, 1, desc );

  REQUIRE( info == 0 );

  // Check correctness

  // A_copy = L**-1 * A_copy
  ptrsm(
    SideFlag::Left, blacspp::Triangle::Lower,
    TransposeFlag::NoTranspose, blacspp::Diagonal::NonUnit,
    M, M, 1., A_SPD_local.data(), 1, 1, desc, 
    A_SPD_copy.data(), 1, 1, desc
  );


  // A_copy = A_copy * L**-H
  ptrsm(
    SideFlag::Right, blacspp::Triangle::Lower,
    TransposeFlag::ConjTranspose, blacspp::Diagonal::NonUnit,
    M, M, 1., A_SPD_local.data(), 1, 1, desc, 
    A_SPD_copy.data(), 1, 1, desc
  );


  std::vector< TestType > gathered;
  if( grid.ipr() == 0 and grid.ipc() == 0 )
    gathered.resize( M*M );

  mat_dist.gather( M, M, gathered.data(), M, A_SPD_copy.data(), M_loc, 0, 0 );

  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    auto tol = M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
    for( auto i = 0; i < M; ++i )
    for( auto j = 0; j < M; ++j ) {
      auto x = std::real( gathered[ i + j*M ] );
      if( i == j ) CHECK( x == Approx( detail::real_t<TestType>(1.) ) );
      else         CHECK( std::abs(x) < tol );
    } 

  }
}
