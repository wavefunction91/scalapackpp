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

  int64_t M = 100, mb = 4;



  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  BlockCyclicMatrix<TestType> A( grid, M, M, mb, mb ),
                              A_SPD( grid, M, M, mb, mb );

  std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );

  for( auto i = 0; i < M; ++i ) 
  if( A.dist().i_own( i, i ) ) {

    auto [ row_idx, col_idx ] = A.dist().local_indx(i, i);
    A.data()[ row_idx + col_idx * A.m_local() ] += 10.;

  }


  // Make A SPD
  pgemm(
    Op::ConjTrans, Op::NoTrans, 
    1., A, A, 0., A_SPD
  );
  auto A_SPD_copy( A_SPD );


  // Perform POTRF
  auto info = ppotrf( blacspp::Uplo::Lower, A_SPD );

  REQUIRE( info == 0 );

  // Check correctness

  // A_copy = L**-1 * A_copy
  ptrsm(
    Side::Left, blacspp::Uplo::Lower,
    Op::NoTrans, blacspp::Diag::NonUnit,
    1., A_SPD, A_SPD_copy
  );


  // A_copy = A_copy * L**-H
  ptrsm(
    Side::Right, blacspp::Uplo::Lower,
    Op::ConjTrans, blacspp::Diag::NonUnit,
    1., A_SPD, A_SPD_copy
  );


  std::vector< TestType > gathered;
  if( grid.ipr() == 0 and grid.ipc() == 0 )
    gathered.resize( M*M );

  A_SPD_copy.gather_from( M, M, gathered.data(), M, 0, 0 );

  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    auto tol = M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
    for( auto i = 0; i < M; ++i )
    for( auto j = 0; j < M; ++j ) {
      auto x = std::real( gathered[ i + j*M ] );
      if( i == j ) CHECK( std::abs(x-TestType(1)) < tol );
      else         CHECK( std::abs(x) < tol );
    } 

  }
}
