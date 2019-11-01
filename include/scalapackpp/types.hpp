#pragma once

#include <blacspp/types.hpp>
#include <array>

namespace scalapackpp {

  using scalapack_int = blacspp::blacs_int;
  using scomplex      = blacspp::scomplex;
  using dcomplex      = blacspp::dcomplex;

  static constexpr scalapack_int scalapack_desc_size = 9;
  using scalapack_desc = std::array< scalapack_int, scalapack_desc_size >;

  enum TransposeFlag {
    NoTranspose,
    Transpose,
    ConjTranspose
  };
}
