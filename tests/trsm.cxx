/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include "ut.hpp"
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/pblas/trsm.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>

#include <blas.hh>

using scalapackpp::internal::scomplex;
using scalapackpp::internal::dcomplex;


SCALAPACKPP_TEST_CASE( "Trsm", "[trsm]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const int64_t M = 100, N = 200;
  int64_t mb = 2, nb = 4;

  std::vector< TestType > A_root, B_ref_root;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A_root.resize( M*M, 1 );
    B_ref_root.resize( M*N, 2 );

    for( auto j = 0; j < M; ++j )
    for( auto i = 0; i < M; ++i )
      if( j > i ) A_root[ i + j*M ] = 0;

    blas::trsm(    
      blas::Layout::ColMajor,
      blas::Side::Left, blas::Uplo::Lower,
      blas::Op::NoTrans, blas::Diag::Unit,
      M, N, 1, A_root.data(), M, B_ref_root.data(), M
    );
  }

  BlockCyclicMatrix<TestType> A( grid, M, M, mb, nb ),
                              B( grid, M, N, mb, nb ),
                              B_ref( grid, M, N, mb, nb );

  std::fill( B.begin(), B.end(), 2. );

  // Scatter triangular matrix and reference
  A.scatter_to( M, M, A_root.data(), M, 0, 0 );
  B_ref.scatter_to( M, N, B_ref_root.data(), M, 0, 0 );


  // Compute with trsm ( B <- A**-1*B )
  ptrsm(
    Side::Left, blacspp::Uplo::Lower,
    Op::NoTrans, blacspp::Diag::Unit,
    1., A, B
  );

  // Check B = C
  for( auto i = 0; i < B.local_size(); ++i )
    CHECK( std::real(B.data()[i]) == Approx( std::real(B_ref.data()[i]) ) );

}

