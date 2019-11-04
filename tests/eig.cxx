#include <catch2/catch.hpp>
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/gemm.hpp>
#include <scalapackpp/syev.hpp>
#include <scalapackpp/heev.hpp>

#include <random>

#define SCALAPACKPP_REAL_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, float, double )

#define SCALAPACKPP_COMPLEX_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, scalapackpp::scomplex, scalapackpp::dcomplex)

#define SCALAPACKPP_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, float, double, scalapackpp::scomplex, scalapackpp::dcomplex)


SCALAPACKPP_REAL_TEST_CASE( "Syev", "[sevp]" ){

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int MB = 4, M = 100;

  auto [M_loc, N_loc] = get_local_dims( grid, M, M, MB, MB, 0, 0);

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 5., 0.25 );

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > Z_local( M_loc * N_loc );
  std::vector< TestType > W( M );

  // Create a symmetric, diagonally dominant matrix
  std::vector< TestType > A;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A.resize( M*M );
    for( auto i = 0; i < M ; ++i )
    for( auto j = 0; j <= i; ++j ) {
      auto x = dist(gen);
      A[ i + j*M ] = x;
      if( i == j ) A[ i + j*M ] += 100;
      else         A[ j + i*M ] = x;
    }
  }

  scatter( grid, M, M, MB, MB, A.data(), M, 0, 0, A_local.data(), M_loc, 0, 0 );

  std::vector< TestType > A_copy( A_local );

  auto context = grid.context();
  auto [desc, i1] = wrappers::descinit( M, M, MB, MB, 0, 0, context, M_loc );
  auto info = psyev(
    VectorFlag::Vectors,
    blacspp::Triangle::Lower,
    M, A_local.data(), 1, 1, desc, W.data(), Z_local.data(), 1, 1, desc
  );

  REQUIRE( info == 0 );

  // Rebuild A
  std::fill( A_local.begin(), A_local.end(), 0 );
  for( auto i = 0; i < M; ++i ) {
    pgemm(
      TransposeFlag::NoTranspose,
      TransposeFlag::Transpose,
      M, M, 1, W[i], Z_local.data(), 1, i+1, desc, Z_local.data(), 1, i+1, desc,
      TestType(1),  A_local.data(), 1, 1, desc
    );
  }

  auto eps = std::numeric_limits<TestType>::epsilon() * 1000;
  for( auto i = 0; i < M_loc*N_loc; ++i ) 
    CHECK( A_local[i] == Approx( A_copy[i] ).epsilon(eps) );
}

SCALAPACKPP_REAL_TEST_CASE( "Syevd", "[sevp]" ){

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int MB = 4, M = 100;

  auto [M_loc, N_loc] = get_local_dims( grid, M, M, MB, MB, 0, 0);

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 5., 0.25 );

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > Z_local( M_loc * N_loc );
  std::vector< TestType > W( M );

  // Create a symmetric, diagonally dominant matrix
  std::vector< TestType > A;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A.resize( M*M );
    for( auto i = 0; i < M ; ++i )
    for( auto j = 0; j <= i; ++j ) {
      auto x = dist(gen);
      A[ i + j*M ] = x;
      if( i == j ) A[ i + j*M ] += 100;
      else         A[ j + i*M ] = x;
    }
  }

  scatter( grid, M, M, MB, MB, A.data(), M, 0, 0, A_local.data(), M_loc, 0, 0 );

  std::vector< TestType > A_copy( A_local );

  auto context = grid.context();
  auto [desc, i1] = wrappers::descinit( M, M, MB, MB, 0, 0, context, M_loc );
  auto info = psyevd(
    VectorFlag::Vectors,
    blacspp::Triangle::Lower,
    M, A_local.data(), 1, 1, desc, W.data(), Z_local.data(), 1, 1, desc
  );

  REQUIRE( info == 0 );

  // Rebuild A
  std::fill( A_local.begin(), A_local.end(), 0 );
  for( auto i = 0; i < M; ++i ) {
    pgemm(
      TransposeFlag::NoTranspose,
      TransposeFlag::Transpose,
      M, M, 1, W[i], Z_local.data(), 1, i+1, desc, Z_local.data(), 1, i+1, desc,
      TestType(1),  A_local.data(), 1, 1, desc
    );
  }

  auto eps = std::numeric_limits<TestType>::epsilon() * 1000;
  for( auto i = 0; i < M_loc*N_loc; ++i ) 
    CHECK( A_local[i] == Approx( A_copy[i] ).epsilon(eps) );
}



SCALAPACKPP_COMPLEX_TEST_CASE( "Heev", "[sevp]" ){

  using namespace scalapackpp;
  using real_t = detail::real_t<TestType>;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int MB = 4, M = 100;

  auto [M_loc, N_loc] = get_local_dims( grid, M, M, MB, MB, 0, 0);

  std::default_random_engine gen;
  std::normal_distribution<real_t> dist( 5., 0.25 );

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > Z_local( M_loc * N_loc );
  std::vector< real_t > W( M );

  // Create a symmetric, diagonally dominant matrix
  std::vector< TestType > A;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A.resize( M*M );
    for( auto i = 0; i < M ; ++i )
    for( auto j = 0; j <= i; ++j ) {
      auto x = dist(gen);
      A[ i + j*M ] = x;
      if( i == j ) A[ i + j*M ] += 100;
      else         A[ j + i*M ] = std::conj(x);
    }
  }

  scatter( grid, M, M, MB, MB, A.data(), M, 0, 0, A_local.data(), M_loc, 0, 0 );

  std::vector< TestType > A_copy( A_local );

  auto context = grid.context();
  auto [desc, i1] = wrappers::descinit( M, M, MB, MB, 0, 0, context, M_loc );
  auto info = pheev(
    VectorFlag::Vectors,
    blacspp::Triangle::Lower,
    M, A_local.data(), 1, 1, desc, W.data(), Z_local.data(), 1, 1, desc
  );

  REQUIRE( info == 0 );

  // Rebuild A
  std::fill( A_local.begin(), A_local.end(), 0 );
  for( auto i = 0; i < M; ++i ) {
    pgemm(
      TransposeFlag::NoTranspose,
      TransposeFlag::ConjTranspose,
      M, M, 1, W[i], Z_local.data(), 1, i+1, desc, Z_local.data(), 1, i+1, desc,
      TestType(1),  A_local.data(), 1, 1, desc
    );
  }

  auto eps = std::numeric_limits<real_t>::epsilon() * 1000;
  for( auto i = 0; i < M_loc*N_loc; ++i ) 
    CHECK( std::real(A_local[i]) == Approx( std::real(A_copy[i]) ).epsilon(eps) );
}

SCALAPACKPP_COMPLEX_TEST_CASE( "Heevd", "[sevp]" ){

  using namespace scalapackpp;
  using real_t = detail::real_t<TestType>;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int MB = 4, M = 100;

  auto [M_loc, N_loc] = get_local_dims( grid, M, M, MB, MB, 0, 0);

  std::default_random_engine gen;
  std::normal_distribution<real_t> dist( 5., 0.25 );

  std::vector< TestType > A_local( M_loc * N_loc );
  std::vector< TestType > Z_local( M_loc * N_loc );
  std::vector< real_t > W( M );

  // Create a symmetric, diagonally dominant matrix
  std::vector< TestType > A;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A.resize( M*M );
    for( auto i = 0; i < M ; ++i )
    for( auto j = 0; j <= i; ++j ) {
      auto x = dist(gen);
      A[ i + j*M ] = x;
      if( i == j ) A[ i + j*M ] += 100;
      else         A[ j + i*M ] = std::conj(x);
    }
  }

  scatter( grid, M, M, MB, MB, A.data(), M, 0, 0, A_local.data(), M_loc, 0, 0 );

  std::vector< TestType > A_copy( A_local );

  auto context = grid.context();
  auto [desc, i1] = wrappers::descinit( M, M, MB, MB, 0, 0, context, M_loc );
  auto info = pheevd(
    VectorFlag::Vectors,
    blacspp::Triangle::Lower,
    M, A_local.data(), 1, 1, desc, W.data(), Z_local.data(), 1, 1, desc
  );

  REQUIRE( info == 0 );

  // Rebuild A
  std::fill( A_local.begin(), A_local.end(), 0 );
  for( auto i = 0; i < M; ++i ) {
    pgemm(
      TransposeFlag::NoTranspose,
      TransposeFlag::ConjTranspose,
      M, M, 1, W[i], Z_local.data(), 1, i+1, desc, Z_local.data(), 1, i+1, desc,
      TestType(1),  A_local.data(), 1, 1, desc
    );
  }

  auto eps = std::numeric_limits<real_t>::epsilon() * 1000;
  for( auto i = 0; i < M_loc*N_loc; ++i ) 
    CHECK( std::real(A_local[i]) == Approx( std::real(A_copy[i]) ).epsilon(eps) );
}
