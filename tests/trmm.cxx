#include <catch2/catch.hpp>
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/trmm.hpp>
#include <scalapackpp/gemm.hpp>
#include <vector>

#define SCALAPACKPP_TEMPLATE_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, float, double, scalapackpp::scomplex, scalapackpp::dcomplex)

SCALAPACKPP_TEMPLATE_TEST_CASE( "Trmm", "[trmm]" ) {

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapackpp::scalapack_int MB = 2, NB = 4;
  const scalapackpp::scalapack_int M = 100, N = 200;

  auto [M_loc1, N_loc1] = scalapackpp::get_local_dims( grid, M, M, MB, MB, 0, 0 );
  auto [M_loc2, N_loc2] = scalapackpp::get_local_dims( grid, M, N, MB, NB, 0, 0 );


  std::vector< TestType > A_root;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A_root.resize( M*M, 1 );
    for( auto j = 0; j < M; ++j )
    for( auto i = 0; i < M; ++i )
      if( j > i ) A_root[ i + j*M ] = 0;
  }

  std::vector< TestType > A_local( M_loc1 * N_loc1 );
  std::vector< TestType > B_local( M_loc2 * N_loc2, TestType(2) );
  std::vector< TestType > C_local( M_loc2 * N_loc2 );

  // Scatter triangular matrix
  scalapackpp::scatter( grid, M, M, MB, MB, A_root.data(), M, 0, 0,
                        A_local.data(), M_loc1, 0, 0 );


  // Compute reference solution with GEMM
  using namespace scalapackpp;
  auto context = grid.context();
  auto [desc_a, i1] = wrappers::descinit( M, M, MB, MB, 0, 0, context, M_loc1 );
  auto [desc_b, i2] = wrappers::descinit( M, N, MB, NB, 0, 0, context, M_loc2 );
  auto [desc_c, i3] = wrappers::descinit( M, N, MB, NB, 0, 0, context, M_loc2 );

  // C <- A*B
  scalapackpp::pgemm( 
    scalapackpp::TransposeFlag::NoTranspose,
    scalapackpp::TransposeFlag::NoTranspose,
    M, N, M, TestType(1), A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b,
    TestType(0), C_local.data(), 1, 1, desc_c 
  );

  // Compute with trmm ( B <- A*B )
  scalapackpp::ptrmm(
    scalapackpp::SideFlag::Left,
    blacspp::Triangle::Lower,
    scalapackpp::TransposeFlag::NoTranspose,
    blacspp::Diagonal::Unit,
    M, N, TestType(1), A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b 
  );

  // Check B = C
  for( auto i = 0; i < C_local.size(); ++i )
    CHECK( std::real(B_local[i]) == Approx( std::real(C_local[i]) ) );

}

