/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/svd.hpp>

SCALAPACKPP_TEST_CASE( "Gesvd", "[svd]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100;
  int64_t N = 10;
  BlockCyclicDist2D mat_dist( grid, 4, 4 );
  int64_t SIZE = std::min(M,N);

  auto [AM_loc, AN_loc]   = mat_dist.get_local_dims( M, N    );
  auto [UM_loc, UN_loc]   = mat_dist.get_local_dims( M, SIZE );
  auto [VTM_loc, VTN_loc] = mat_dist.get_local_dims( SIZE, N );

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A_local( AM_loc * AN_loc );
  std::vector< TestType > U_local( UM_loc * UN_loc );
  std::vector< TestType > VT_local( VTM_loc * VTN_loc );

  std::vector< detail::real_t<TestType> > S( SIZE );

  std::generate( A_local.begin(), A_local.end(), [&](){ return dist(gen); } );
  auto A_copy = A_local;

  auto desc_a  = mat_dist.descinit_noerror( M, N,    AM_loc  );
  auto desc_u  = mat_dist.descinit_noerror( M, SIZE, UM_loc  );
  auto desc_vt = mat_dist.descinit_noerror( SIZE, N, VTM_loc );

  pgesvd( VectorFlag::Vectors, VectorFlag::Vectors, M, N,
    A_local.data(), 1, 1, desc_a, S.data(),
    U_local.data(), 1, 1, desc_u, VT_local.data(), 1, 1, desc_vt );

  // Rebuild A
  std::fill( A_local.begin(), A_local.end(), 0 );
  for( auto i = 0; i < SIZE; ++i ) {
    pgemm(
      TransposeFlag::NoTranspose,
      TransposeFlag::NoTranspose,
      M, N, 1, S[i], 
      U_local.data(),  1, i+1, desc_u, 
      VT_local.data(), i+1, 1, desc_vt,
      TestType(1.),  A_local.data(), 1, 1, desc_a
    );
  }

  auto eps = std::numeric_limits<detail::real_t<TestType>>::epsilon() * M*M*N;
  for( auto i = 0; i < AM_loc*AN_loc; ++i ) 
    CHECK( std::real(A_local[i]) == Approx( std::real(A_copy[i]) ).epsilon(eps) );
}
