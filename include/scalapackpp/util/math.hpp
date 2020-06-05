/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <cstdlib>
#include <type_traits>

namespace scalapackpp::detail {

template <typename M, typename N>
std::make_signed_t<std::common_type_t<std::remove_reference_t<M>,std::remove_reference_t<N>>> div_ceil( M m, N n ) {

  using ct = std::make_signed_t<std::common_type_t<M,N>>;

  ct m_c = m;
  ct n_c = n;

  auto d = std::div( m, n );
  return d.quot + !!d.rem;
  
}

}
