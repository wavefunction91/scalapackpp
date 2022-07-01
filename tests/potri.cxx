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







SCALAPACKPP_REAL_TEST_CASE( "Potri", "[potri]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100, mb = 4;

  std::default_random_engine gen(mpi.rank());
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

  info = ppotri( blacspp::Uplo::Lower, A_SPD );

  REQUIRE( info == 0 );

  // Symmetrize
  std::copy( A_SPD.begin(), A_SPD.end(), A.begin() );

  // Zero out upper triangle
  fill_triangle( blacspp::Uplo::Upper, A,     0. ); 
  fill_triangle( blacspp::Uplo::Upper, A_SPD, 0. ); 

  pgeadd( Op::ConjTrans, 1., A, 1., A_SPD );

  // Account for double counting of diagonals
  for( auto i = 0; i < M; ++i )
  if( A.dist().i_own( i, i )  ) {
    auto [ I, J ] = A.dist().local_indx( i, i );
    A_SPD.data()[ I + J*A.m_local()] /= 2.;
  }


  // Check correctness

  pgemm(
    Op::NoTrans, Op::NoTrans, 
    1., A_SPD, A_SPD_copy, 0., A
  );


  auto tol = 100*M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
  const detail::real_t<TestType> one = 1.;
  for( auto i = 0; i < M; ++i )
  for( auto j = 0; j < M; ++j )
  if( A.dist().i_own( i, j ) ) {

    auto [ I, J ] = A.dist().local_indx( i, j );
    auto x = std::real( A.data()[ I + J * A.m_local() ] );

    if( i == j ) CHECK( x == Approx( one ).epsilon(tol) );
    else         CHECK( std::abs(x) < tol );

  }

}
