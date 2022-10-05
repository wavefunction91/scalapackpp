/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/eigenvalue_problem/sevp.hpp>

template <typename T>
T smart_conj( T x ) {
  if constexpr (std::is_floating_point_v<T>)
    return x;
  else return std::conj(x);
}




SCALAPACKPP_TEST_CASE( "SEVP Drivers", "[sevp]" ){

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100;
  int64_t mb = 4;

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 5., 0.25 );

  BlockCyclicMatrix<TestType> A_sca( grid, M, M, mb, mb ),
                              Z_sca( grid, M, M, mb, mb );

  std::vector< detail::real_t<TestType> > W( M );

  // Create a symmetric, diagonally dominant matrix
  std::vector< TestType > A;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A.resize( M*M );
    for( auto i = 0; i < M ; ++i )
    for( auto j = 0; j <= i; ++j ) {
      auto x = dist(gen);
      A[ i + j*M ] = x;
      if( i == j ) A[ i + j*M ] += 100;
      else         A[ j + i*M ] = smart_conj(x);
    }
  }

  A_sca.scatter_to( M, M, A.data(), M, 0, 0 );

  auto A_copy = A_sca;

  SECTION( "HEREIG" ) {
    auto info = hereig(
      Job::Vec, blacspp::Uplo::Lower,
      A_sca, W.data(), Z_sca
    );
    REQUIRE( info == 0 );
  }
  SECTION( "HEREIGD" ) {
    auto info = hereigd(
      Job::Vec, blacspp::Uplo::Lower,
      A_sca, W.data(), Z_sca
    );
    REQUIRE( info == 0 );
  }

  // Rebuild A
  std::fill( A_sca.begin(), A_sca.end(), 0 );
  for( auto i = 0; i < M; ++i ) {
    pgemm(
      Op::NoTrans,
      Op::ConjTrans,
      M, M, 1, W[i], Z_sca.data(), 1, i+1, Z_sca.desc(), Z_sca.data(), 1, i+1, Z_sca.desc(),
      1,  A_sca.data(), 1, 1, A_sca.desc()
    );
  }

  auto eps = std::numeric_limits<detail::real_t<TestType>>::epsilon() * 1000;
  for( auto i = 0; i < A_sca.local_size(); ++i ) 
    CHECK( std::real(A_sca.data()[i]) == Approx( std::real(A_copy.data()[i]) ).epsilon(eps) );
  MPI_Barrier(MPI_COMM_WORLD);
}
