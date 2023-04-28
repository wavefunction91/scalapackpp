/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/pblas/symv.hpp>
#include <scalapackpp/util/type_conversions.hpp>

#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T>
  psymv( Uplo uplo, 
         int64_t N, detail::type_identity_t<T> ALPHA, 
         const T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         const T* X, int64_t IX, int64_t JX, const scalapack_desc& DESCX, int64_t INCX,
         detail::type_identity_t<T> BETA,
         T* Y, int64_t IY, int64_t JY, const scalapack_desc& DESCY, int64_t INCY ) {

  assert( X != Y );

  auto UPLO = char( uplo );

  wrappers::psymv( &UPLO, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX, INCX,
                   BETA, Y, IY, JY, DESCY, INCY ); 

}

}
