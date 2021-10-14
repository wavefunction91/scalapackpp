/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/lacpy.hpp>
#include <scalapackpp/util/type_conversions.hpp>

#include <scalapackpp/block_cyclic_matrix.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  placpy( Uplo uplo, int64_t M, int64_t N, const T* A, int64_t IA, int64_t JA,
          const scalapack_desc& DESCA, T* B, int64_t IB, int64_t JB,
          const scalapack_desc& DESCB ) {

  auto UPLO = char( uplo );
  wrappers::placpy( &UPLO, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  placpy( Uplo uplo, const BlockCyclicMatrix<T>& A, BlockCyclicMatrix<T>&B ) {

  placpy( uplo, A.m(), A.n(), A.data(), 1, 1, A.desc(), B.data(), 1, 1, B.desc() );

}

}
