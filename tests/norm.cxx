/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/matrix_norm/general_norm.hpp>
#include <scalapackpp/matrix_norm/symmetric_norm.hpp>

#include <lapack.hh>



SCALAPACKPP_TEST_CASE( "GeneralNorm", "[norm]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100;

  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );
  auto desc = mat_dist.descinit_noerror( M, M, M_loc );

  std::default_random_engine gen;
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

  std::vector< TestType > A_local( M_loc * N_loc );
  std::generate( A_local.begin(), A_local.end(), [&](){ return dist(gen); } );



  std::vector< TestType > A_gathered;
  if( mpi.rank() == 0 ) A_gathered.resize( M*M );

  mat_dist.gather(M, M, A_gathered.data(), M, A_local.data(), M_loc, 0, 0);

  detail::real_t<TestType> f_ref, i_ref, o_ref, m_ref;
  if( mpi.rank() == 0 ) {
    f_ref = lapack::lange( lapack::Norm::Fro, M, M, A_gathered.data(), M );
    i_ref = lapack::lange( lapack::Norm::Inf, M, M, A_gathered.data(), M );
    o_ref = lapack::lange( lapack::Norm::One, M, M, A_gathered.data(), M );
    m_ref = std::abs(std::real(*std::max_element( A_gathered.begin(), A_gathered.end(),
              [](const auto a, const auto b){ return std::abs(a) < std::abs(b); }
            )));

  }

  auto dtype =
    std::is_same_v<detail::real_t<TestType>,float> ? MPI_FLOAT : MPI_DOUBLE;

  MPI_Bcast( &f_ref, 1, dtype, 0, MPI_COMM_WORLD );
  MPI_Bcast( &i_ref, 1, dtype, 0, MPI_COMM_WORLD );
  MPI_Bcast( &o_ref, 1, dtype, 0, MPI_COMM_WORLD );
  MPI_Bcast( &m_ref, 1, dtype, 0, MPI_COMM_WORLD );

  auto f_calc = general_norm( mat_dist, FrobeniusNorm, M, M, A_local.data(), 1, 1,
    desc );
  auto i_calc = general_norm( mat_dist, InfinityNorm, M, M, A_local.data(), 1, 1,
    desc );
  auto o_calc = general_norm( mat_dist, OneNorm, M, M, A_local.data(), 1, 1,
    desc );
  auto m_calc = general_norm( mat_dist, AbsMax, M, M, A_local.data(), 1, 1,
    desc );


  CHECK( f_calc == Approx( f_ref ) );
  CHECK( i_calc == Approx( i_ref ) );
  CHECK( o_calc == Approx( o_ref ) );
  CHECK( m_calc == Approx( m_ref ) );

  MPI_Barrier(MPI_COMM_WORLD);
}



#if 0
SCALAPACKPP_TEST_CASE( "SymmetricNorm", "[norm]" ) {

  using namespace scalapackpp;
  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100;


  std::vector< TestType > A_gathered;


  detail::real_t<TestType> f_ref, i_ref, o_ref, m_ref;
  if( mpi.rank() == 0 ) {
    A_gathered.resize(M*M);

    std::default_random_engine gen;
    std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );

    for( int j = 0; j < M; ++j )
    for( int i = j; i < M; ++i ) {
      auto v = dist(gen);
      A_gathered[ i + j*M ] = v;
      A_gathered[ j + i*M ] = v;
    }

    f_ref = lange( 'F', M, M, A_gathered.data(), M );
    i_ref = lange( 'I', M, M, A_gathered.data(), M );
    o_ref = lange( 'O', M, M, A_gathered.data(), M );
    m_ref = std::abs(std::real(*std::max_element( A_gathered.begin(), A_gathered.end(),
              [](const auto a, const auto b){ return std::abs(a) < std::abs(b); }
            )));

  }

  auto dtype =
    std::is_same_v<detail::real_t<TestType>,float> ? MPI_FLOAT : MPI_DOUBLE;

  MPI_Bcast( &f_ref, 1, dtype, 0, MPI_COMM_WORLD );
  MPI_Bcast( &i_ref, 1, dtype, 0, MPI_COMM_WORLD );
  MPI_Bcast( &o_ref, 1, dtype, 0, MPI_COMM_WORLD );
  MPI_Bcast( &m_ref, 1, dtype, 0, MPI_COMM_WORLD );

  BlockCyclicDist2D mat_dist( grid, 4, 4 );

  auto [M_loc, N_loc] = mat_dist.get_local_dims( M, M );
  auto desc = mat_dist.descinit_noerror( M, M, M_loc );

  std::vector< TestType > A_local( M_loc * N_loc );
  mat_dist.scatter(M, M, A_gathered.data(), M, A_local.data(), M_loc, 0, 0);




  auto f_calc = symmetric_norm( mat_dist, FrobeniusNorm, blacspp::Triangle::Lower ,M, A_local.data(), 1, 1,
    desc );                                              
 // auto i_calc = symmetric_norm( mat_dist, InfinityNorm,  blacspp::Triangle::Lower ,M, A_local.data(), 1, 1,
 //   desc );                                              
 // auto o_calc = symmetric_norm( mat_dist, OneNorm,       blacspp::Triangle::Lower ,M, A_local.data(), 1, 1,
 //   desc );                                              
  auto m_calc = symmetric_norm( mat_dist, AbsMax,        blacspp::Triangle::Lower ,M, A_local.data(), 1, 1,
    desc );


  CHECK( f_calc == Approx( f_ref ) );
 // CHECK( i_calc == Approx( i_ref ) );
 // CHECK( o_calc == Approx( o_ref ) );
  CHECK( m_calc == Approx( m_ref ) );

}
#endif
