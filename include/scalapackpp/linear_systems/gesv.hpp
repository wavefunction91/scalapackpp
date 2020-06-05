/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/linear_systems/gesv.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  pgesv( scalapack_int N, scalapack_int NRHS, 
    T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
    scalapack_int* IPIV,
    T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB ) {

  return wrappers::pgesv( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB );

}
          
}
