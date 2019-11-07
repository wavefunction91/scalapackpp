#pragma once
#include <scalapackpp/wrappers/geadd.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T, typename ALPHAT, typename BETAT>
std::enable_if_t<
  detail::scalapack_supported_v<T> and 
  std::is_convertible_v<ALPHAT,T> and
  std::is_convertible_v<BETAT,T>
>
  pgeadd( TransposeFlag trans, scalapack_int M, scalapack_int N, ALPHAT ALPHA,
        const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
        BETAT BETA, 
        T* C, scalapack_int IC, scalapack_int JC, const scalapack_desc& DESCC ) {

  auto TRANS = detail::type_string( trans );

  const T ALPHA_t = T(ALPHA);
  const T BETA_t  = T(BETA);  


  wrappers::pgeadd( TRANS.c_str(), M, N, ALPHA_t, A, IA, JA, DESCA, BETA_t,
                    C, IC, JC, DESCC );

}

}
