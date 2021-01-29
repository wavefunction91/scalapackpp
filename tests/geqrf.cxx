/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/factorizations/geqrf.hpp>
#include <scalapackpp/factorizations/geqpf.hpp>
#include <scalapackpp/householder/generate_q_householder.hpp>
#include <scalapackpp/pblas/gemm.hpp>

SCALAPACKPP_TEST_CASE( "Geqrf", "[geqrf]" ) {

  using namespace scalapackpp;

  std::shared_ptr<const blacspp::Grid> grid = std::make_shared<const blacspp::Grid>(blacspp::Grid::square_grid( MPI_COMM_WORLD ));
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100, N = 50;
  int64_t mb = 4;

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A;

  if( grid->ipr() == 0 and grid->ipc() == 0 ) {
    A.resize( M*N );
    std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );
  }

  BlockCyclicMatrix<TestType> A_sca( grid, M, N, mb, mb );
  std::vector< TestType > 
    TAU_local( local_col_from_desc( std::min(M,N), A_sca.desc() ) );

  A_sca.scatter_to( M, N, A.data(), M, 0, 0 );
  auto A_cpy(A_sca);

  auto info = pgeqrf( A_sca, TAU_local.data() );
  REQUIRE( info == 0 );

  // Generate Q
  auto Q_sca(A_sca);
  info = generate_q_householder( N, Q_sca, TAU_local.data() );
  REQUIRE( info == 0 );

  // Compute R = Q**H * A
  BlockCyclicMatrix<TestType> R_sca( grid, N, N, mb, mb );
  pgemm( Op::ConjTrans, Op::NoTrans,
         1., Q_sca, A_cpy, 0., R_sca );



  // Check correctnesss
  std::vector<TestType> R;
  if( grid->ipc() == 0 and grid->ipr() == 0 ) {
    R.resize( N*N );
  }

  R_sca.gather_from( N, N, R.data(), N, 0, 0 );
  A_sca.gather_from( M, N, A.data(), M, 0, 0 );

  if( grid->ipc() == 0 and grid->ipr() == 0 ) {

    const auto tol = M*N*std::numeric_limits<detail::real_t<TestType>>::epsilon();
    for( int j = 0; j <  N; ++j )
    for( int i = 0; i <= j; ++i ) {
      CHECK( std::abs(R[i + j*N] - A[ i + j*M ]) < tol );
    }

  }


}


SCALAPACKPP_TEST_CASE( "Geqpf", "[geqpf]" ) {

  using namespace scalapackpp;

  std::shared_ptr<const blacspp::Grid> grid = std::make_shared<const blacspp::Grid>(blacspp::Grid::square_grid( MPI_COMM_WORLD ));
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100, N = 50;
  int64_t mb = 4;

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A;

  if( grid->ipr() == 0 and grid->ipc() == 0 ) {
    A.resize( M*N );
    std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );
  }

  BlockCyclicMatrix<TestType> A_sca( grid, M, N, mb, mb );
  std::vector< TestType > 
    TAU_local( local_col_from_desc( std::min(M,N), A_sca.desc() ) );
  std::vector< int64_t  > IPIV_local( local_col_from_desc( N, A_sca.desc() ) );

  A_sca.scatter_to( M, N, A.data(), M, 0, 0 );
  auto A_cpy(A_sca);

  auto info = pgeqpf( A_sca, IPIV_local.data(), TAU_local.data() );
  REQUIRE( info == 0 );

  // TODO: check for correctness
}
