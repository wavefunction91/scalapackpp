/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/pblas/trmm.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#include <blacspp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T, typename ALPHAT>
detail::enable_if_t<
  detail::scalapack_supported<T>::value and 
  std::is_convertible<ALPHAT,T>::value
>
  ptrmm( Side side, blacspp::Triangle uplo, Op trans, 
         blacspp::Diagonal diag,
         int64_t M, int64_t N, ALPHAT ALPHA, 
         const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {

  auto SIDE = char( side );
  auto UPLO = blacspp::detail::type_string( uplo );
  auto TRANS = char( trans );
  auto DIAG  = blacspp::detail::type_string( diag );

  const T ALPHA_t = T(ALPHA);

  wrappers::ptrmm( &SIDE, UPLO.c_str(), &TRANS, DIAG.c_str(),
                   M, N, ALPHA_t, A, IA, JA, DESCA, B, IB, JB, DESCB );

}

}
