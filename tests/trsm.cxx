/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/trsm.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>

#include <blas.hh>

using scalapackpp::scalapack_int;
using scalapackpp::scomplex;
using scalapackpp::dcomplex;


SCALAPACKPP_TEST_CASE( "Trsm", "[trsm]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapack_int M = 100, N = 200;

  BlockCyclicDist2D mat_dist( grid, 2, 4 );

  auto [M_loc1, N_loc1] = mat_dist.get_local_dims( M, M );
  auto [M_loc2, N_loc2] = mat_dist.get_local_dims( M, N );


  std::vector< TestType > A_root, B_ref_root;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A_root.resize( M*M, 1 );
    B_ref_root.resize( M*N, 2 );

    for( auto j = 0; j < M; ++j )
    for( auto i = 0; i < M; ++i )
      if( j > i ) A_root[ i + j*M ] = 0;

    blas::trsm(    
      blas::Layout::ColMajor,
      blas::Side::Left, blas::Uplo::Lower,
      blas::Op::NoTrans, blas::Diag::Unit,
      M, N, 1, A_root.data(), M, B_ref_root.data(), M
    );
  }

  std::vector< TestType > A_local( M_loc1*N_loc1 );
  std::vector< TestType > B_local( M_loc2*N_loc2, 2 );
  std::vector< TestType > B_ref_local( M_loc2*N_loc2 );

  // Scatter triangular matrix and reference
  mat_dist.scatter( M, M, A_root.data(), M, A_local.data(), M_loc1, 0, 0 );
  mat_dist.scatter( M, N, B_ref_root.data(), M, B_ref_local.data(), M_loc2, 0, 0 );


  // Compute reference solution with GEMM
  auto desc_a = mat_dist.descinit_noerror( M, M, M_loc1 );
  auto desc_b = mat_dist.descinit_noerror( M, N, M_loc2 );

  // Compute with trsm ( B <- A**-1*B )
  ptrsm(
    SideFlag::Left, blacspp::Triangle::Lower,
    TransposeFlag::NoTranspose, blacspp::Diagonal::Unit,
    M, N, 1, A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b 
  );

  // Check B = C
  for( auto i = 0; i < B_local.size(); ++i )
    CHECK( std::real(B_local[i]) == Approx( std::real(B_ref_local[i]) ) );

}

