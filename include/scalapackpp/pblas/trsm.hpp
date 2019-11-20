#pragma once
#include <scalapackpp/wrappers/pblas/trsm.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T, typename ALPHAT>
std::enable_if_t<
  detail::scalapack_supported_v<T> and 
  std::is_convertible_v<ALPHAT,T>
>
  ptrsm( SideFlag side, blacspp::Triangle uplo, TransposeFlag trans, 
         blacspp::Diagonal diag,
         scalapack_int M, scalapack_int N, ALPHAT ALPHA, 
         const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  auto SIDE = detail::type_string( side );
  auto UPLO = blacspp::detail::type_string( uplo );
  auto TRANS = detail::type_string( trans );
  auto DIAG  = blacspp::detail::type_string( diag );

  const T ALPHA_t = T(ALPHA);

  wrappers::ptrsm( SIDE.c_str(), UPLO.c_str(), TRANS.c_str(), DIAG.c_str(),
                   M, N, ALPHA_t, A, IA, JA, DESCA, B, IB, JB, DESCB );

}

}
