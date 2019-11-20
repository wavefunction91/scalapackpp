#pragma once
#include <scalapackpp/wrappers/linear_systems/trtrs.hpp>
#include <blacspp/util/type_conversions.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  ptrtrs( blacspp::Triangle uplo, TransposeFlag trans, blacspp::Diagonal diag,
    scalapack_int N, scalapack_int NRHS, 
    const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  auto UPLO = blacspp::detail::type_string( uplo );
  auto DIAG = blacspp::detail::type_string( diag );
  auto TRANS = scalapackpp::detail::type_string( trans );
  return wrappers::ptrtrs( UPLO.c_str(), TRANS.c_str(), DIAG.c_str(), N, NRHS,
           A, IA, JA, DESCA, B, IB, JB, DESCB );

}
          
}
