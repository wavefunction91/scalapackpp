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

  const int64_t M = 100, N = 200;
  const int64_t mb = 4;



  BlockCyclicMatrix<TestType> A( grid, M, M, mb, mb ),
                              B( grid, M, N, mb, mb );

  std::fill( A.begin(), A.end(), 1 );
  std::fill( B.begin(), B.end(), 2 );
  auto B_ref( B );


  // Zero out upper triangle
  fill_triangle( blacspp::Uplo::Upper, A, 0. ); 


  // Compute reference sol with trsm
  ptrsm(
    Side::Left, blacspp::Uplo::Lower, Op::NoTrans,
    blacspp::Diag::Unit, 1., A, B_ref );

  // Compute with trtrs ( B <- A**-1*B )
  auto info = ptrtrs(
    blacspp::Uplo::Lower,
    Op::NoTrans, blacspp::Diag::Unit,
    A, B
  );

  REQUIRE( info == 0 );

  // Check B = C
  for( auto i = 0; i < B.local_size(); ++i )
    CHECK( std::real(B.data()[i]) == Approx( std::real(B_ref.data()[i]) ) );

  MPI_Barrier(MPI_COMM_WORLD);
}

