/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/svd.hpp>

SCALAPACKPP_TEST_CASE( "Gesvd", "[svd]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  int64_t M = 100;
  int64_t N = 10;
  int64_t SIZE = std::min(M,N);
  int64_t mb = 4;

  BlockCyclicMatrix<TestType> A( grid, M, N, mb, mb ),
                              U( grid, M, SIZE, mb, mb ),
                              VT( grid, SIZE, N, mb, mb );
  std::vector< detail::real_t<TestType> > S( SIZE );


  std::default_random_engine gen(mpi.rank());
  std::normal_distribution<detail::real_t<TestType>> dist( 0., 1. );
  std::generate( A.begin(), A.end(), [&](){ return dist(gen); } );
  auto A_copy = A;


  pgesvd( Job::Vec, Job::Vec, A, S.data(), U, VT );


  // Rebuild A
  std::fill( A.begin(), A.end(), 0 );
  for( auto i = 0; i < SIZE; ++i ) {
    pgemm(
      Op::NoTrans,
      Op::NoTrans,
      M, N, 1, S[i], 
      U.data(),  1, i+1, U.desc(), 
      VT.data(), i+1, 1, VT.desc(),
      1.,  A.data(), 1, 1, A.desc()
    );
  }

  auto eps = std::numeric_limits<detail::real_t<TestType>>::epsilon() * M*M*N;
  for( auto i = 0; i < A.local_size(); ++i ) 
    CHECK( std::real(A.data()[i]) == Approx( std::real(A_copy.data()[i]) ).epsilon(eps) );
}
