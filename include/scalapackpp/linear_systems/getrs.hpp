#pragma once
#include <scalapackpp/wrappers/linear_systems/getrs.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  pgetrs( TransposeFlag trans, scalapack_int N, scalapack_int NRHS, 
    const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    scalapack_int* IPIV,
    T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  auto TRANS = detail::type_string( trans );
  return wrappers::pgetrs( TRANS.c_str(), N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB );

}
          
}
