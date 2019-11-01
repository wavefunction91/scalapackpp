#include <catch2/catch.hpp>
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/gemm.hpp>
#include <vector>

#define SCALAPACKPP_TEMPLATE_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, float, double, scalapackpp::scomplex, scalapackpp::dcomplex)

SCALAPACKPP_TEMPLATE_TEST_CASE( "Gemm", "[gemm]" ) {

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapackpp::scalapack_int MB = 4, NB = 4;
  const scalapackpp::scalapack_int M = MB * grid.npr(), N = NB * grid.npr();


  std::vector< TestType > A_local( MB*MB, TestType(2) );
  std::vector< TestType > B_local( MB*NB, TestType(3) );
  std::vector< TestType > C_local( MB*NB, TestType(5) );

  using namespace scalapackpp;

  auto context = grid.context();
  auto [desc_a, i1] = wrappers::descinit( M, M, MB, MB, 0, 0, context, MB );
  auto [desc_b, i2] = wrappers::descinit( M, N, MB, NB, 0, 0, context, MB );
  auto [desc_c, i3] = wrappers::descinit( M, N, MB, NB, 0, 0, context, MB );


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
