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
#include <scalapackpp/linear_systems/potrs.hpp>







SCALAPACKPP_TEST_CASE( "Potrs", "[potrs]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100, NRHS = 20;
  int64_t mb = 4;


  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  BlockCyclicMatrix<TestType> A( grid, M, M,    mb, mb ),
                              B( grid, M, NRHS, mb, mb ),
                              A_SPD( grid, M, M, mb, mb );

  std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );
  std::generate( B.begin(), B.end(), [&](){ return dist(gen); } );
  auto B_copy( B );

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

  // Solve linear system
  auto info = ppotrf( blacspp::Uplo::Lower, A_SPD );
  REQUIRE( info == 0 );
  info = ppotrs( blacspp::Uplo::Lower, A_SPD, B );
  REQUIRE( info == 0 );

  // Check correctness
  pgemm( Op::NoTrans, Op::NoTrans, 
         -1., A_SPD_copy, B, 1., B_copy );

  auto tol = 1000*M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
  for( auto x : B_copy ) CHECK( std::abs(x) < tol );

  MPI_Barrier(MPI_COMM_WORLD);
}
