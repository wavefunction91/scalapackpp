#include <catch2/catch.hpp>
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/descinit.hpp>
#include <scalapackpp/trsm.hpp>
#include <scalapackpp/gemm.hpp>
#include <scalapackpp/potrf.hpp>

#include <random>


#define SCALAPACKPP_TEMPLATE_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, float, double, scalapackpp::scomplex, scalapackpp::dcomplex)




SCALAPACKPP_TEMPLATE_TEST_CASE( "Potrf", "[potrf]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapackpp::scalapack_int MB = 4, M = 100;

  auto [M_loc, N_loc] = scalapackpp::get_local_dims( grid, M, M, MB, MB, 0, 0 );

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > A_SPD_local( M_loc * N_loc );

  std::generate( A_local.begin(), A_local.end(), [&](){ return dist(gen); } );

  for( auto i = 0; i < M; ++i ) {
    if( grid.ipr() == ( i / MB ) % grid.npr() and
        grid.ipc() == ( i / MB ) % grid.npc() ) {

      auto block_row = i / (MB * grid.npr());
      auto block_col = i / (MB * grid.npc());

      auto row_idx = block_row * MB + i % MB;
      auto col_idx = block_col * MB + i % MB;


      A_local[ row_idx + col_idx * M_loc ] += 10.;
    }
  }

  auto desc = descinit_noerror( grid, M, M, MB, MB, 0, 0, M_loc );

  // Make A SPD
  scalapackpp::pgemm(
    scalapackpp::TransposeFlag::ConjTranspose,
    scalapackpp::TransposeFlag::NoTranspose,
    M, M, M, TestType(1.), A_local.data(), 1, 1, desc, A_local.data(), 1, 1, desc,
    TestType(0.), A_SPD_local.data(), 1, 1, desc
  );
  std::vector< TestType > A_SPD_copy( A_SPD_local );


  // Perform POTRF
  auto info = ppotrf( blacspp::Triangle::Lower, M, A_SPD_local.data(), 1, 1, desc );

  REQUIRE( info == 0 );

  // Check correctness

  // A_copy = L**-1 * A_copy
  scalapackpp::ptrsm(
    scalapackpp::SideFlag::Left,
    blacspp::Triangle::Lower,
    scalapackpp::TransposeFlag::NoTranspose,
    blacspp::Diagonal::NonUnit,
    M, M, TestType(1.), A_SPD_local.data(), 1, 1, desc, 
    A_SPD_copy.data(), 1, 1, desc
  );


  // A_copy = A_copy * L**-H
  scalapackpp::ptrsm(
    scalapackpp::SideFlag::Right,
    blacspp::Triangle::Lower,
    scalapackpp::TransposeFlag::ConjTranspose,
    blacspp::Diagonal::NonUnit,
    M, M, TestType(1.), A_SPD_local.data(), 1, 1, desc, 
    A_SPD_copy.data(), 1, 1, desc
  );


  std::vector< TestType > gathered;
  if( grid.ipr() == 0 and grid.ipc() == 0 )
    gathered.resize( M*M );

  scalapackpp::gather( grid, M, M, MB, MB, gathered.data(), M, 0, 0,
                       A_SPD_copy.data(), M_loc, 0, 0 );

  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    for( auto i = 0; i < M; ++i )
    for( auto j = 0; j < M; ++j ) {
      auto x = std::real( gathered[ i + j*M ] );
      if( i == j ) CHECK( x == Approx( detail::real_t<TestType>(1.) ) );
      else         CHECK( std::abs(x) < M*M*std::numeric_limits<detail::real_t<TestType>>::epsilon() );
    } 

  }
}
