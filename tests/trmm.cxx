#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/trmm.hpp>
#include <scalapackpp/gemm.hpp>


SCALAPACKPP_TEST_CASE( "Trmm", "[trmm]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapack_int M = 100, N = 200;

  BlockCyclicDist2D mat_dist( grid, 2, 4 );

  auto [M_loc1, N_loc1] = mat_dist.get_local_dims( M, M );
  auto [M_loc2, N_loc2] = mat_dist.get_local_dims( M, N );


  std::vector< TestType > A_root;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A_root.resize( M*M, 1 );
    for( auto j = 0; j < M; ++j )
    for( auto i = 0; i < M; ++i )
      if( j > i ) A_root[ i + j*M ] = 0;
  }

  std::vector< TestType > A_local( M_loc1 * N_loc1 );
  std::vector< TestType > B_local( M_loc2 * N_loc2, 2 );
  std::vector< TestType > C_local( M_loc2 * N_loc2 );

  // Scatter triangular matrix
  mat_dist.scatter( M, M, A_root.data(), M, A_local.data(), M_loc1, 0, 0 );


  // Compute reference solution with GEMM
  auto desc_a = mat_dist.descinit_noerror( M, M, M_loc1 );
  auto desc_b = mat_dist.descinit_noerror( M, N, M_loc2 );

  // C <- A*B
  pgemm( 
    TransposeFlag::NoTranspose, TransposeFlag::NoTranspose, M, N, M, 
    1, A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b,
    0, C_local.data(), 1, 1, desc_b
  );

  // Compute with trmm ( B <- A*B )
  ptrmm(
    SideFlag::Left, blacspp::Triangle::Lower,
    TransposeFlag::NoTranspose, blacspp::Diagonal::Unit,
    M, N, 1, A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b 
  );

  // Check B = C
  for( auto i = 0; i < C_local.size(); ++i )
    CHECK( std::real(B_local[i]) == Approx( std::real(C_local[i]) ) );

}

