/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/geadd.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T, typename ALPHAT, typename BETAT>
std::enable_if_t<
  detail::scalapack_supported<T>::value and 
  std::is_convertible<ALPHAT,T>::value and
  std::is_convertible<BETAT,T>::value
>
  pgeadd( TransposeFlag trans, int64_t M, int64_t N, ALPHAT ALPHA,
        const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
        BETAT BETA, 
        T* C, int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {

  assert( A != C );

  auto TRANS = detail::type_string( trans );

  const T ALPHA_t = T(ALPHA);
  const T BETA_t  = T(BETA);  


  wrappers::pgeadd( TRANS.c_str(), M, N, ALPHA_t, A, IA, JA, DESCA, BETA_t,
                    C, IC, JC, DESCC );

}

}
