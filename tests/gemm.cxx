#include <catch2/catch.hpp>
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/gemm.hpp>
#include <vector>

#define SCALAPACKPP_TEMPLATE_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, float, double, scalapackpp::scomplex, scalapackpp::dcomplex)

SCALAPACKPP_TEMPLATE_TEST_CASE( "Gemm", "[gemm]" ) {

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapackpp::scalapack_int MB = 2, NB = 4;
  const scalapackpp::scalapack_int M = 100, N = 200;

  auto [M_loc1, N_loc1] = scalapackpp::get_local_dims( grid, M, M, MB, MB, 0, 0 );
  auto [M_loc2, N_loc2] = scalapackpp::get_local_dims( grid, M, N, MB, NB, 0, 0 );


  std::vector< TestType > A_local( M_loc1 * N_loc1, TestType(2) );
  std::vector< TestType > B_local( M_loc2 * N_loc2, TestType(3) );
  std::vector< TestType > C_local( M_loc2 * N_loc2, TestType(5) );

  using namespace scalapackpp;

  auto context = grid.context();
  auto [desc_a, i1] = wrappers::descinit( M, M, MB, MB, 0, 0, context, M_loc1 );
  auto [desc_b, i2] = wrappers::descinit( M, N, MB, NB, 0, 0, context, M_loc2 );
  auto [desc_c, i3] = wrappers::descinit( M, N, MB, NB, 0, 0, context, M_loc2 );


  scalapackpp::pgemm( 
    scalapackpp::TransposeFlag::NoTranspose,
    scalapackpp::TransposeFlag::NoTranspose,
    M, N, M, TestType(1.), A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b,
    TestType(2.), C_local.data(), 1, 1, desc_c 
  );


  const TestType val = TestType(2 * 5 + M*2*3);
  for( auto x : C_local )
    CHECK( std::real(x) == Approx( std::real(val) ) );

}
