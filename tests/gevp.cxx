/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/eigenvalue_problem/gevp.hpp>

template <typename T>
T smart_conj( T x ) {
  if constexpr (std::is_floating_point_v<T>)
    return x;
  else return std::conj(x);
}


SCALAPACKPP_TEST_CASE( "Hereig_gen", "[gevp]" ){

  using namespace scalapackpp;

  using real_t = detail::real_t<TestType>;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int M = 100;
  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );

  std::default_random_engine gen;
  std::normal_distribution<real_t> dist( 5., 0.25 );

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > B_local( M_loc * N_loc );
  std::vector< TestType > Z_local( M_loc * N_loc );
  std::vector< real_t > W( M );

  // Create a symmetric, diagonally dominant matrix (A)
  // and SPD matrix B
  std::vector< TestType > A, B;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A.resize( M*M );
    B.resize( M*M );
    for( auto i = 0; i < M ; ++i )
    for( auto j = 0; j <= i; ++j ) {
      auto x = dist(gen);
      A[ i + j*M ] = x;
      if( i == j ) A[ i + j*M ] += 100;
      else         A[ j + i*M ] = smart_conj(x);

      auto y = dist(gen);
      B[ i + j*M ] = y;
      if( i == j ) B[ i + j*M ] += 10;
      else         B[ j + i*M ] = smart_conj(y);
    }
  }



  mat_dist.scatter( M, M, A.data(), M, A_local.data(), M_loc, 0, 0 );
  mat_dist.scatter( M, M, B.data(), M, Z_local.data(), M_loc, 0, 0 );

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );

  // Make B SPD
  pgemm(
    TransposeFlag::ConjTranspose, TransposeFlag::NoTranspose, M, M, M,
    1., Z_local.data(), 1, 1, desc, Z_local.data(), 1, 1, desc,
    0., B_local.data(), 1, 1, desc 
  );

  std::vector< TestType > A_copy( A_local );
  std::vector< TestType > B_copy( B_local );
 

  auto info = hereig_gen(
    VectorFlag::Vectors,
    blacspp::Triangle::Lower,
    M, A_local.data(), 1, 1, desc, B_local.data(), 1, 1, desc, 
    W.data(), Z_local.data(), 1, 1, desc
  );

  REQUIRE( info == 0 );

  // Rebuild A
  std::fill( A_local.begin(), A_local.end(), 0 );

  // B_local <- B * Z
  pgemm( TransposeFlag::NoTranspose, TransposeFlag::NoTranspose,
         M, M, M, 1., B_copy.data(), 1, 1, desc, Z_local.data(), 1, 1, desc,
         0., B_local.data(), 1, 1, desc );

  // Z = B_local
  Z_local = std::move(B_local);

  // A <- X * E * X**H
  for( auto i = 0; i < M; ++i ) {
    pgemm(
      TransposeFlag::NoTranspose,
      TransposeFlag::ConjTranspose,
      M, M, 1, W[i], Z_local.data(), 1, i+1, desc, Z_local.data(), 1, i+1, desc,
      TestType(1),  A_local.data(), 1, 1, desc
    );
  }




  auto eps = std::numeric_limits<real_t>::epsilon() * 10000;
  for( auto i = 0; i < M_loc*N_loc; ++i ) 
    CHECK( std::real(A_local[i]) == Approx( std::real(A_copy[i]) ).epsilon(eps) );
}

