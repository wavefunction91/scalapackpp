/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic_matrix.hpp>
#include <scalapackpp/pblas/inner_product.hpp>

SCALAPACKPP_REAL_TEST_CASE("Dot", "[inner_product]") {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t N = 200;
  const int64_t nb = 4;

  BlockCyclicMatrix<TestType> 
    V(grid, N, N, nb, nb), W(grid, N, N, nb, nb );

  std::vector<TestType> V_loc(N*N), W_loc(N*N);
  std::iota(V_loc.begin(), V_loc.end(), 0);
  std::iota(W_loc.begin(), W_loc.end(), 1);

  V.scatter_to( N, N, V_loc.data(), N, 0, 0 );
  W.scatter_to( N, N, W_loc.data(), N, 0, 0 );

  for(int i = 0; i < N; ++i)
  for(int j = 0; j < N; ++j) {

    TestType ref_inner = 0.0;
    for( int k = 0; k < N; ++k )
      ref_inner += V_loc[k + i*N] * W_loc[k + j*N];

    auto inner = pdot( N, V.data(), 1, i+1, V.desc(), 1, W.data(), 1, j+1,
      W.desc(), 1 ); 

    auto own_col_i = (i / nb) % grid.npc();
    auto own_col_j = (j / nb) % grid.npc();

    if(grid.ipc() == own_col_i or grid.ipc() == own_col_j)
      REQUIRE(inner == Approx(ref_inner));

    grid.barrier(blacspp::Scope::All);
  }
  
}

SCALAPACKPP_COMPLEX_TEST_CASE("Dotu", "[inner_product]") {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t N = 200;
  const int64_t nb = 4;

  BlockCyclicMatrix<TestType> 
    V(grid, N, N, nb, nb), W(grid, N, N, nb, nb );

  std::vector<TestType> V_loc(N*N), W_loc(N*N);
  std::iota(V_loc.begin(), V_loc.end(), 0);
  std::iota(W_loc.begin(), W_loc.end(), 1);
  for( auto i = 0; i < N*N; ++i ) {
    auto p = TestType(std::real(V_loc[i]),  std::real(W_loc[i]));
    auto m = TestType(std::real(W_loc[i]), -std::real(V_loc[i]));
    V_loc[i] = p;
    W_loc[i] = m;
  }

  V.scatter_to( N, N, V_loc.data(), N, 0, 0 );
  W.scatter_to( N, N, W_loc.data(), N, 0, 0 );

  for(int i = 0; i < N; ++i)
  for(int j = 0; j < N; ++j) {

    TestType ref_inner = 0.0;
    for( int k = 0; k < N; ++k )
      ref_inner += V_loc[k + i*N] * W_loc[k + j*N];

    auto inner = pdotu( N, V.data(), 1, i+1, V.desc(), 1, W.data(), 1, j+1,
      W.desc(), 1 ); 

    auto own_col_i = (i / nb) % grid.npc();
    auto own_col_j = (j / nb) % grid.npc();

    if(grid.ipc() == own_col_i or grid.ipc() == own_col_j) {
      REQUIRE(std::real(inner) == Approx(std::real(ref_inner)));
      REQUIRE(std::imag(inner) == Approx(std::imag(ref_inner)));
    }

    grid.barrier(blacspp::Scope::All);
  }
  
}


SCALAPACKPP_COMPLEX_TEST_CASE("Dotc", "[inner_product]") {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t N = 200;
  const int64_t nb = 4;

  BlockCyclicMatrix<TestType> 
    V(grid, N, N, nb, nb), W(grid, N, N, nb, nb );

  std::vector<TestType> V_loc(N*N), W_loc(N*N);
  std::iota(V_loc.begin(), V_loc.end(), 0);
  std::iota(W_loc.begin(), W_loc.end(), 1);
  for( auto i = 0; i < N*N; ++i ) {
    auto p = TestType(std::real(V_loc[i]),  std::real(W_loc[i]));
    auto m = TestType(std::real(W_loc[i]), -std::real(V_loc[i]));
    V_loc[i] = p;
    W_loc[i] = m;
  }

  V.scatter_to( N, N, V_loc.data(), N, 0, 0 );
  W.scatter_to( N, N, W_loc.data(), N, 0, 0 );

  for(int i = 0; i < N; ++i)
  for(int j = 0; j < N; ++j) {

    TestType ref_inner = 0.0;
    for( int k = 0; k < N; ++k )
      ref_inner += std::conj(V_loc[k + i*N]) * W_loc[k + j*N];

    auto inner = pdotc( N, V.data(), 1, i+1, V.desc(), 1, W.data(), 1, j+1,
      W.desc(), 1 ); 

    auto own_col_i = (i / nb) % grid.npc();
    auto own_col_j = (j / nb) % grid.npc();

    if(grid.ipc() == own_col_i or grid.ipc() == own_col_j) {
      REQUIRE(std::real(inner) == Approx(std::real(ref_inner)));
      REQUIRE(std::imag(inner) == Approx(std::imag(ref_inner)));
    }

    grid.barrier(blacspp::Scope::All);
  }
  
}
