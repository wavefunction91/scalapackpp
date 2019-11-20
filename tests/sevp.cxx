#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/eigenvalue_problem/sevp.hpp>

template <typename T>
T smart_conj( T x ) {
  if constexpr (std::is_floating_point_v<T>)
    return x;
  else return std::conj(x);
}



SCALAPACKPP_REAL_TEST_CASE( "Syev", "[sevp]" ){

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int M = 100;
  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );

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

  mat_dist.scatter( M, M, A.data(), M, A_local.data(), M_loc, 0, 0 );

  std::vector< TestType > A_copy( A_local );

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );
  auto info = psyev(
    VectorFlag::Vectors, blacspp::Triangle::Lower,
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

  scalapack_int M = 100;
  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );

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

  mat_dist.scatter( M, M, A.data(), M, A_local.data(), M_loc, 0, 0 );

  std::vector< TestType > A_copy( A_local );

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );
  auto info = psyevd(
    VectorFlag::Vectors, blacspp::Triangle::Lower,
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

  scalapack_int M = 100;
  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );

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
      else         A[ j + i*M ] = smart_conj(x);
    }
  }

  mat_dist.scatter( M, M, A.data(), M, A_local.data(), M_loc, 0, 0 );

  std::vector< TestType > A_copy( A_local );

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );
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

  scalapack_int M = 100;
  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );

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
      else         A[ j + i*M ] = smart_conj(x);
    }
  }

  mat_dist.scatter( M, M, A.data(), M, A_local.data(), M_loc, 0, 0 );

  std::vector< TestType > A_copy( A_local );

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );
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















SCALAPACKPP_TEST_CASE( "Hereig", "[sevp]" ){

  using namespace scalapackpp;

  using real_t = detail::real_t<TestType>;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int M = 100;
  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );

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
      else         A[ j + i*M ] = smart_conj(x);
    }
  }

  mat_dist.scatter( M, M, A.data(), M, A_local.data(), M_loc, 0, 0 );

  std::vector< TestType > A_copy( A_local );

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );
  auto info = hereig(
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

SCALAPACKPP_TEST_CASE( "Hereigd", "[sevp]" ){

  using namespace scalapackpp;

  using real_t = detail::real_t<TestType>;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  scalapack_int M = 100;
  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );

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
      else         A[ j + i*M ] = smart_conj(x);
    }
  }

  mat_dist.scatter( M, M, A.data(), M, A_local.data(), M_loc, 0, 0 );

  std::vector< TestType > A_copy( A_local );

  auto desc = mat_dist.descinit_noerror( M, M, M_loc );
  auto info = hereigd(
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
