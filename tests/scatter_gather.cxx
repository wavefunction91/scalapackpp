#include <catch2/catch.hpp>
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/information.hpp>
#include <blacspp/information.hpp>
#include <vector>
#include <iostream>

#define SCALAPACKPP_TEMPLATE_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, float, double, scalapackpp::scomplex, scalapackpp::dcomplex)


SCALAPACKPP_TEMPLATE_TEST_CASE( "Gather", "[scatter-gather]" ) {

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapackpp::scalapack_int MB = 2, NB = 4;
  const scalapackpp::scalapack_int M = 100, N = 200;

  auto [M_loc, N_loc] = scalapackpp::get_local_dims( grid, M, N, MB, NB, 0, 0 );

  std::vector< TestType > data_local( M_loc * N_loc, TestType(mpi.rank()) );
  std::vector< TestType > data_gathered;
  if( grid.ipr() == 0 and grid.ipc() == 0 )
    data_gathered.resize( M*N, TestType(-1) );


  scalapackpp::gather( grid, M, N, MB, NB, data_gathered.data(), M, 
                       0, 0, data_local.data(), M_loc, 0, 0 );

  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    for( int i = 0; i < M; ++i )
    for( int j = 0; j < N; ++j ) {
      auto row = (i / MB) % grid.npr();
      auto col = (j / NB) % grid.npc();

      auto rank = blacspp::coordinate_rank( grid, row, col );
      CHECK( data_gathered[i + j*M] == TestType(rank) );
    }

  } else {
    CHECK( data_gathered.size() == 0 );
  }

  for( auto x : data_local ) CHECK( x == TestType(mpi.rank()) );

  grid.barrier( blacspp::Scope::All );

}



SCALAPACKPP_TEMPLATE_TEST_CASE( "Scatter", "[scatter-gather]" ) {

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapackpp::scalapack_int MB = 2, NB = 4;
  const scalapackpp::scalapack_int M = 100, N = 200;

  auto [M_loc, N_loc] = scalapackpp::get_local_dims( grid, M, N, MB, NB, 0, 0 );


  std::vector< TestType > data_local( M_loc * N_loc, TestType(-1) );
  std::vector< TestType > data_gathered;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {

    data_gathered.resize( M*N );

    for( int i = 0; i < M; ++i )
    for( int j = 0; j < N; ++j ) {
      auto row = (i / MB) % grid.npr();
      auto col = (j / NB) % grid.npc();
      data_gathered[i + j*M] = blacspp::coordinate_rank( grid, row, col );
    }

  }


  scalapackpp::scatter( grid, M, N, MB, NB, data_gathered.data(), M, 
                        0, 0, data_local.data(), M_loc, 0, 0 );

  for( auto x : data_local ) CHECK( x == TestType( mpi.rank() ) );

}
