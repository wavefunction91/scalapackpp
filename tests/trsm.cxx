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

using scalapackpp::scalapack_int;
using scalapackpp::scomplex;
using scalapackpp::dcomplex;

extern "C" {

void strsm_( const char*, const char*, const char*, const char*, 
             const scalapack_int*, const scalapack_int*, const float*,
             const float*, const scalapack_int*, float*, const scalapack_int* );
void dtrsm_( const char*, const char*, const char*, const char*, 
             const scalapack_int*, const scalapack_int*, const double*,
             const double*, const scalapack_int*, double*, const scalapack_int* );
void ctrsm_( const char*, const char*, const char*, const char*, 
             const scalapack_int*, const scalapack_int*, const scomplex*,
             const scomplex*, const scalapack_int*, scomplex*, const scalapack_int* );
void ztrsm_( const char*, const char*, const char*, const char*, 
             const scalapack_int*, const scalapack_int*, const dcomplex*,
             const dcomplex*, const scalapack_int*, dcomplex*, const scalapack_int* );
}

namespace scalapackpp {

template <typename T, typename ALPHAT>
void trsm( SideFlag side, blacspp::Triangle uplo, TransposeFlag trans, blacspp::Diagonal diag,
         scalapack_int M, scalapack_int N, ALPHAT ALPHA, 
         const T* A, scalapack_int LDA, T* B, scalapack_int LDB ) {

  auto SIDE = detail::type_string( side );
  auto UPLO = blacspp::detail::type_string( uplo );
  auto TRANS = detail::type_string( trans );
  auto DIAG  = blacspp::detail::type_string( diag );

  T ALPHA_t = ALPHA;

  if constexpr ( std::is_same_v< T, float > )
  strsm_( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
          &M, &N, &ALPHA_t, A, &LDA, B, &LDB );
  if constexpr ( std::is_same_v< T, double > )
  dtrsm_( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
          &M, &N, &ALPHA_t, A, &LDA, B, &LDB );
  if constexpr ( std::is_same_v< T, scomplex > )
  ctrsm_( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
          &M, &N, &ALPHA_t, A, &LDA, B, &LDB );
  if constexpr ( std::is_same_v< T, dcomplex > )
  ztrsm_( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
          &M, &N, &ALPHA_t, A, &LDA, B, &LDB );

}

}

SCALAPACKPP_TEST_CASE( "Trsm", "[trsm]" ) {

  using namespace scalapackpp;

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapack_int M = 100, N = 200;

  BlockCyclicDist2D mat_dist( grid, 2, 4 );

  auto [M_loc1, N_loc1] = mat_dist.get_local_dims( M, M );
  auto [M_loc2, N_loc2] = mat_dist.get_local_dims( M, N );


  std::vector< TestType > A_root, B_ref_root;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A_root.resize( M*M, 1 );
    B_ref_root.resize( M*N, 2 );

    for( auto j = 0; j < M; ++j )
    for( auto i = 0; i < M; ++i )
      if( j > i ) A_root[ i + j*M ] = 0;

    scalapackpp::trsm(    
      SideFlag::Left, blacspp::Triangle::Lower,
      TransposeFlag::NoTranspose, blacspp::Diagonal::Unit,
      M, N, 1, A_root.data(), M, B_ref_root.data(), M
    );
  }

  std::vector< TestType > A_local( M_loc1*N_loc1 );
  std::vector< TestType > B_local( M_loc2*N_loc2, 2 );
  std::vector< TestType > B_ref_local( M_loc2*N_loc2 );

  // Scatter triangular matrix and reference
  mat_dist.scatter( M, M, A_root.data(), M, A_local.data(), M_loc1, 0, 0 );
  mat_dist.scatter( M, N, B_ref_root.data(), M, B_ref_local.data(), M_loc2, 0, 0 );


  // Compute reference solution with GEMM
  auto desc_a = mat_dist.descinit_noerror( M, M, M_loc1 );
  auto desc_b = mat_dist.descinit_noerror( M, N, M_loc2 );

  // Compute with trsm ( B <- A**-1*B )
  ptrsm(
    SideFlag::Left, blacspp::Triangle::Lower,
    TransposeFlag::NoTranspose, blacspp::Diagonal::Unit,
    M, N, 1, A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b 
  );

  // Check B = C
  for( auto i = 0; i < B_local.size(); ++i )
    CHECK( std::real(B_local[i]) == Approx( std::real(B_ref_local[i]) ) );

}

