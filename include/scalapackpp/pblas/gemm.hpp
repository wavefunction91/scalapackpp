#pragma once
#include <scalapackpp/wrappers/pblas/gemm.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T, typename ALPHAT, typename BETAT>
std::enable_if_t<
  detail::scalapack_supported_v<T> and 
  std::is_convertible_v<ALPHAT,T> and
  std::is_convertible_v<BETAT,T>
>
  pgemm( TransposeFlag transa, TransposeFlag transb,
         scalapack_int M, scalapack_int N, scalapack_int K, ALPHAT ALPHA, 
         const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         const T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
         BETAT BETA,
         T* C, scalapack_int IC, scalapack_int JC, const scalapack_desc& DESCC ) {

  assert( A != C );
  assert( B != C );

  auto TRANSA = detail::type_string( transa );
  auto TRANSB = detail::type_string( transb );

  const T ALPHA_t = T(ALPHA);
  const T BETA_t  = T(BETA);  

  wrappers::pgemm( TRANSA.c_str(), TRANSB.c_str(), M, N, K, ALPHA_t, A, IA, JA,
                   DESCA, B, IB, JB, DESCB, BETA_t, C, IC, JC, DESCC ); 

}


}
