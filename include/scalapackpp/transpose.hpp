#pragma once
#include <scalapackpp/geadd.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  transpose( scalapack_int M, scalapack_int N, 
    const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB 
  ) {

  pgeadd( TransposeFlag::Transpose, M, N, 1., A, IA, JA, DESCA, 
          0., B, IB, JB, DESCB );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  conj_transpose( scalapack_int M, scalapack_int N, 
    const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB 
  ) {

  if constexpr ( detail::scalapack_real_supported_v<T> )
    transpose( M, N, A, IA, JA, DESCA, B, IB, JB, DESCB );
  else
    pgeadd( TransposeFlag::ConjTranspose, M, N, 1., A, IA, JA, DESCA, 
            0., B, IB, JB, DESCB );

}
          

} 
