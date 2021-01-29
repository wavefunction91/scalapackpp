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
#include <scalapackpp/matrix_inverse/getri.hpp>


SCALAPACKPP_TEST_CASE( "Getri", "[getri]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100, mb = 4;


  std::default_random_engine gen(mpi.rank());
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  BlockCyclicMatrix<TestType> A( grid, M, M, mb, mb );

  std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );

  for( auto i = 0; i < M; ++i ) 
  if( A.dist().i_own( i, i ) ) {

    auto [ row_idx, col_idx ] = A.dist().local_indx(i, i);
    A.data()[ row_idx + col_idx * A.m_local() ] += 10.;

  }

  auto A_inv( A );
  std::vector<int64_t> IPIV( A.m_local() + A.dist().mb() );
  auto info = pgetrf( A_inv, IPIV.data() );
  REQUIRE( info == 0 );

  info = pgetri( A_inv, IPIV.data() );
  REQUIRE( info == 0 );


  // Check correctness
  BlockCyclicMatrix<TestType> I_approx( grid, M, M, mb, mb );
  pgemm(
    TransposeFlag::NoTranspose, TransposeFlag::NoTranspose, 
    1., A, A_inv, 0., I_approx
  );


  auto tol = 10*M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
  const detail::real_t<TestType> one = 1.;
  for( auto i = 0; i < M; ++i )
  for( auto j = 0; j < M; ++j )
  if( A.dist().i_own( i, j ) ) {

    auto [ I, J ] = A.dist().local_indx( i, j );
    auto x = std::real( I_approx.data()[ I + J * A.m_local() ] );

    if( i == j ) CHECK( x == Approx( one ).epsilon(tol) );
    else         CHECK( std::abs(x) < tol );

  }

}
