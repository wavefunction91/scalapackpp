#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/matrix_inverse/trtri.hpp>
#include <scalapackpp/fill_triangle.hpp>

#if 0
SCALAPACKPP_REAL_TEST_CASE( "trtri", "[trtri]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int M = 100;

  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A_local( M_loc * N_loc );

  std::generate( A_local.begin(), A_local.end(), [&](){ return dist(gen); } );

  for( auto i = 0; i < M; ++i ) 
  if( mat_dist.i_own( i, i ) ) {

    auto [ row_idx, col_idx ] = mat_dist.local_indx(i, i);
    A_local[ row_idx + col_idx * M_loc ] += 10.;

  }

  // Zero upper triangle
  fill_triangle( mat_dist, blacspp::Triangle::Upper, M, M,
    A_local.data(), M_loc, 0. );

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );

  std::vector< TestType > A_inv_local( A_local );
  auto info = ptrtri( blacspp::Triangle::Lower, blacspp::Diagonal::NonUnit, M,
    A_local.data(), 1, 1, desc );
  REQUIRE( info == 0 );


  // Check correctness
  decltype(A_local) I_approx( M_loc * N_loc );
  pgemm(
    TransposeFlag::NoTranspose, TransposeFlag::NoTranspose, M, M, M, 
    1., A_local.data(), 1, 1, desc, A_inv_local.data(), 1, 1, desc,
    0., I_approx.data(), 1, 1, desc
  );


  auto tol = 10*M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon();
  const detail::real_t<TestType> one = 1.;
  for( auto i = 0; i < M; ++i )
  for( auto j = 0; j < M; ++j )
  if( mat_dist.i_own( i, j ) ) {

    auto [ I, J ] = mat_dist.local_indx( i, j );
    auto x = std::real( I_approx[ I + J * M_loc ] );

    if( i == j ) CHECK( x == Approx( one ).epsilon(tol) );
    else         CHECK( std::abs(x) < tol );

  }

}
#endif
