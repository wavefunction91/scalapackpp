#pragma once
#include <scalapackpp/wrappers/pblas/scal.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T, typename ALPHAT>
std::enable_if_t<
  detail::scalapack_supported_v<T> and 
  std::is_convertible_v<ALPHAT,T>
>
  pscal( scalapack_int N, ALPHAT ALPHA, T* X, scalapack_int IX, scalapack_int JX, 
         const scalapack_desc& DESCX, scalapack_int INCX ) {

  const T ALPHA_t = T(ALPHA);
  wrappers::pscal( N, ALPHA_t, X, IX, JX, DESCX, INCX );

}


}

