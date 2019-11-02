#pragma once
#include <scalapackpp/wrappers/trmm.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  ptrmm( SideFlag side, blacspp::Triangle uplo, TransposeFlag trans, blacspp::Diagonal diag,
         scalapack_int M, scalapack_int N, T ALPHA, 
         const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  auto SIDE = detail::type_string( side );
  auto UPLO = blacspp::detail::type_string( uplo );
  auto TRANS = detail::type_string( trans );
  auto DIAG  = blacspp::detail::type_string( diag );

  wrappers::ptrmm( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
                   M, N, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB );

}

}
