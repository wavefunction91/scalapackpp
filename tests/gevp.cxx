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


SCALAPACKPP_TEST_CASE( "GEVP Drivers", "[gevp]" ){

  using namespace scalapackpp;

  using real_t = detail::real_t<TestType>;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100;
  int64_t mb = 4;

  std::default_random_engine gen;
  std::normal_distribution<real_t> dist( 5., 0.25 );

  BlockCyclicMatrix<TestType> A_sca( grid, M, M, mb, mb ),
                              B_sca( grid, M, M, mb, mb ),
                              Z_sca( grid, M, M, mb, mb );

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



  A_sca.scatter_to( M, M, A.data(), M, 0, 0 );
  Z_sca.scatter_to( M, M, B.data(), M, 0, 0 );

  // Make B SPD
  pgemm(
    Op::ConjTrans, Op::NoTrans, 
    1., Z_sca, Z_sca, 0., B_sca
  );

  auto A_copy( A_sca );
  auto B_copy( B_sca );
 
  SECTION( "HEREIG_GEN" ) {

    blacspp::Triangle uplo;
    SECTION( "Lower" ) { uplo = blacspp::Triangle::Lower; }
    SECTION( "Upper" ) { uplo = blacspp::Triangle::Upper; }

    auto info = hereig_gen(
      Job::Vec,
      uplo,
      A_sca, B_sca, W.data(), Z_sca
    );

    REQUIRE( info == 0 );

  }

  SECTION( "HEREIGD_GEN" ) {

    blacspp::Triangle uplo;
    SECTION( "Lower" ) { uplo = blacspp::Triangle::Lower; }
    SECTION( "Upper" ) { uplo = blacspp::Triangle::Upper; }

    auto info = hereigd_gen(
      Job::Vec,
      uplo,
      A_sca, B_sca, W.data(), Z_sca
    );

    REQUIRE( info == 0 );

  }

  // Rebuild A
  std::fill( A_sca.begin(), A_sca.end(), 0 );

  // B_local <- B * Z
  pgemm( Op::NoTrans, Op::NoTrans,
         1., B_copy, Z_sca, 0., B_sca );


  // A <- X * E * X**H
  for( auto i = 0; i < M; ++i ) {
    pgemm(
      Op::NoTrans,
      Op::ConjTrans,
      M, M, 1, W[i], B_sca.data(), 1, i+1, B_sca.desc(), B_sca.data(), 1, i+1, B_sca.desc(),
      1,  A_sca.data(), 1, 1, A_sca.desc()
    );
  }




  auto eps = std::numeric_limits<real_t>::epsilon() * 10000;
  for( auto i = 0; i < A_sca.local_size(); ++i ) 
    CHECK( std::real(A_sca.data()[i]) == Approx( std::real(A_copy.data()[i]) ).epsilon(eps) );
}

