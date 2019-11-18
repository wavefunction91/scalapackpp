#pragma once
#include <scalapackpp/wrappers/potrs.hpp>
#include <blacspp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  ppotrs( blacspp::Triangle uplo, scalapack_int N, scalapack_int NRHS, 
    const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {


  auto UPLO = blacspp::detail::type_string( uplo );
  return wrappers::ppotrs( UPLO.c_str(), N, NRHS, A, IA, JA, DESCA,
                           B, IB, JB, DESCB );
  
}

}
