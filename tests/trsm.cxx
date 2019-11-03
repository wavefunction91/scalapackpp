#include <catch2/catch.hpp>
#include <scalapackpp/scatter_gather.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/trsm.hpp>
#include <vector>

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

template <typename T>
void trsm( SideFlag side, blacspp::Triangle uplo, TransposeFlag trans, blacspp::Diagonal diag,
         scalapack_int M, scalapack_int N, T ALPHA, 
         const T* A, scalapack_int LDA, T* B, scalapack_int LDB ) {

  auto SIDE = detail::type_string( side );
  auto UPLO = blacspp::detail::type_string( uplo );
  auto TRANS = detail::type_string( trans );
  auto DIAG  = blacspp::detail::type_string( diag );

  if constexpr ( std::is_same_v< T, float > )
  strsm_( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
          &M, &N, &ALPHA, A, &LDA, B, &LDB );
  if constexpr ( std::is_same_v< T, double > )
  dtrsm_( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
          &M, &N, &ALPHA, A, &LDA, B, &LDB );
  if constexpr ( std::is_same_v< T, scomplex > )
  ctrsm_( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
          &M, &N, &ALPHA, A, &LDA, B, &LDB );
  if constexpr ( std::is_same_v< T, dcomplex > )
  ztrsm_( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
          &M, &N, &ALPHA, A, &LDA, B, &LDB );

}

}

#define SCALAPACKPP_TEMPLATE_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, float, double, scalapackpp::scomplex, scalapackpp::dcomplex)

SCALAPACKPP_TEMPLATE_TEST_CASE( "Trsm", "[trsm]" ) {

  blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
  blacspp::mpi_info mpi( MPI_COMM_WORLD );

  const scalapackpp::scalapack_int MB = 2, NB = 4;
  const scalapackpp::scalapack_int M = 100, N = 200;

  auto [M_loc1, N_loc1] = scalapackpp::get_local_dims( grid, M, M, MB, MB, 0, 0 );
  auto [M_loc2, N_loc2] = scalapackpp::get_local_dims( grid, M, N, MB, NB, 0, 0 );


  std::vector< TestType > A_root, B_ref_root;
  if( grid.ipr() == 0 and grid.ipc() == 0 ) {
    A_root.resize( M*M, 1 );
    B_ref_root.resize( M*N, 2 );

    for( auto j = 0; j < M; ++j )
    for( auto i = 0; i < M; ++i )
      if( j > i ) A_root[ i + j*M ] = 0;

    scalapackpp::trsm(    
      scalapackpp::SideFlag::Left,
      blacspp::Triangle::Lower,
      scalapackpp::TransposeFlag::NoTranspose,
      blacspp::Diagonal::Unit,
      M, N, TestType(1), A_root.data(), M, B_ref_root.data(), M
    );
  }

  std::vector< TestType > A_local( M_loc1*N_loc1 );
  std::vector< TestType > B_local( M_loc2*N_loc2, TestType(2) );
  std::vector< TestType > B_ref_local( M_loc2*N_loc2 );

  // Scatter triangular matrix and reference
  scalapackpp::scatter( grid, M, M, MB, MB, A_root.data(), M, 0, 0,
                        A_local.data(), M_loc1, 0, 0 );
  scalapackpp::scatter( grid, M, N, MB, NB, B_ref_root.data(), M, 0, 0,
                        B_ref_local.data(), M_loc2, 0, 0 );


  // Compute reference solution with GEMM
  using namespace scalapackpp;
  auto context = grid.context();
  auto [desc_a, i1] = wrappers::descinit( M, M, MB, MB, 0, 0, context, M_loc1 );
  auto [desc_b, i2] = wrappers::descinit( M, N, MB, NB, 0, 0, context, M_loc2 );

  // Compute with trsm ( B <- A**-1*B )
  scalapackpp::ptrsm(
    scalapackpp::SideFlag::Left,
    blacspp::Triangle::Lower,
    scalapackpp::TransposeFlag::NoTranspose,
    blacspp::Diagonal::Unit,
    M, N, TestType(1), A_local.data(), 1, 1, desc_a, B_local.data(), 1, 1, desc_b 
  );

  // Check B = C
  for( auto i = 0; i < B_local.size(); ++i )
    CHECK( std::real(B_local[i]) == Approx( std::real(B_ref_local[i]) ) );

}

