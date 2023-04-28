/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/wrappers/pblas/nrm2.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T,detail::real_t<T>>
  pnrm2( int64_t N, const T* X, int64_t IX, int64_t JX, const scalapack_desc& DESCX, 
        int64_t INCX ) { 

  detail::real_t<T> nrm;
  wrappers::pnrm2( N, &nrm, X, IX, JX, DESCX, INCX );
  return nrm;

}

}

