/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/factorizations/potrf.hpp>
#include <scalapackpp/matrix_inverse/potri.hpp>
#include <scalapackpp/geadd.hpp>
#include <scalapackpp/fill_triangle.hpp>







SCALAPACKPP_TEST_CASE( "Potri", "[potri]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100;

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

  info = ppotri( blacspp::Triangle::Lower, M, A_SPD_local.data(), 1, 1, desc );

  REQUIRE( info == 0 );

  // Symmetrize
  std::copy( A_SPD_local.begin(), A_SPD_local.end(), A_local.begin() );

  // Zero out upper triangle
  fill_triangle( mat_dist, blacspp::Triangle::Upper, M, M, 
    A_local.data(), M_loc, 0. );
  fill_triangle( mat_dist, blacspp::Triangle::Upper, M, M, 
    A_SPD_local.data(), M_loc, 0. );

  pgeadd( TransposeFlag::ConjTranspose, M, M, 1., A_local.data(), 1, 1, desc,
          1., A_SPD_local.data(), 1, 1, desc );

  // Account for double counting of diagonals
  for( auto i = 0; i < M; ++i )
  if( mat_dist.i_own( i, i )  ) {
    auto [ I, J ] = mat_dist.local_indx( i, i );
    A_SPD_local[ I + J*M_loc] /= 2.;
  }


  // Check correctness

  pgemm(
    TransposeFlag::NoTranspose, TransposeFlag::NoTranspose, M, M, M, 
    1., A_SPD_local.data(), 1, 1, desc, A_SPD_copy.data(), 1, 1, desc,
    0., A_local.data(), 1, 1, desc
  );


  auto tol = 10*M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
  const detail::real_t<TestType> one = 1.;
  for( auto i = 0; i < M; ++i )
  for( auto j = 0; j < M; ++j )
  if( mat_dist.i_own( i, j ) ) {

    auto [ I, J ] = mat_dist.local_indx( i, j );
    auto x = std::real( A_local[ I + J * M_loc ] );

    if( i == j ) CHECK( x == Approx( one ).epsilon(tol) );
    else         CHECK( std::abs(x) < tol );

  }

}
