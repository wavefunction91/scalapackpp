#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/linear_systems/posv.hpp>







SCALAPACKPP_TEST_CASE( "Posv", "[posv]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int M = 100, NRHS = 20;

  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );
  auto [xxxxx, NRHS_loc ] = mat_dist.get_local_dims( M, NRHS );

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > B_local( M_loc * NRHS_loc );
  std::vector< TestType > A_SPD_local( M_loc * N_loc );

  std::generate( A_local.begin(), A_local.end(), [&](){ return dist(gen); } );
  std::generate( B_local.begin(), B_local.end(), [&](){ return dist(gen); } );
  auto B_copy( B_local );

  for( auto i = 0; i < M; ++i ) 
  if( mat_dist.i_own( i, i ) ) {

    auto [ row_idx, col_idx ] = mat_dist.local_indx(i, i);
    A_local[ row_idx + col_idx * M_loc ] += 10.;

  }

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );
  auto desc_rhs = mat_dist.descinit_noerror( M, NRHS, M_loc );

  // Make A SPD
  pgemm(
    TransposeFlag::ConjTranspose, TransposeFlag::NoTranspose, M, M, M, 
    1., A_local.data(), 1, 1, desc, A_local.data(), 1, 1, desc,
    0., A_SPD_local.data(), 1, 1, desc
  );
  std::vector< TestType > A_SPD_copy( A_SPD_local );

  // Solve linear system
  auto info = pposv( blacspp::Triangle::Lower, M, NRHS, A_SPD_local.data(), 1, 1, desc,
                 B_local.data(), 1, 1, desc_rhs );

  REQUIRE( info == 0 );
  // Check correctness
  pgemm( TransposeFlag::NoTranspose, TransposeFlag::NoTranspose, M, NRHS, M,
         -1., A_SPD_copy.data(), 1, 1, desc, B_local.data(), 1, 1, desc_rhs,
         1. , B_copy.data(), 1, 1, desc_rhs );

  auto tol = 100*M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
  for( auto x : B_copy ) CHECK( std::abs(x) < tol );

}
