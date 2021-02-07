/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/pblas/symm.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T, typename ALPHAT, typename BETAT>
detail::enable_if_t<
  detail::scalapack_real_supported<T>::value and 
  std::is_convertible<ALPHAT,T>::value and
  std::is_convertible<BETAT,T>::value
>
  psymm( Side side, Uplo uplo,
         int64_t M, int64_t N, ALPHAT ALPHA, 
         const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         const T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB,
         BETAT BETA,
         T* C, int64_t IC, int64_t JC, const scalapack_desc& DESCC ) {

  assert( A != C );
  assert( B != C );

  auto SIDE = char( side );
  auto UPLO = char( uplo );

  const T ALPHA_t = T(ALPHA);
  const T BETA_t  = T(BETA);  

  wrappers::psymm( &SIDE, &UPLO, M, N, ALPHA_t, A, IA, JA,
                   DESCA, B, IB, JB, DESCB, BETA_t, C, IC, JC, DESCC ); 

}


}
