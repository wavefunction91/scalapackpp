#pragma once

#include <blacspp/types.hpp>
#include <array>

namespace scalapackpp {

  using scalapack_int = blacspp::blacs_int;
  using scomplex      = blacspp::scomplex;
  using dcomplex      = blacspp::dcomplex;

  static constexpr scalapack_int scalapack_desc_size = 9;
  using scalapack_desc = std::array< scalapack_int, scalapack_desc_size >;


  struct block_cyclic_coordinate {
    scalapack_int process_id;
    scalapack_int block_idx;
    scalapack_int local_idx;
  };


  enum TransposeFlag {
    NoTranspose,
    Transpose,
    ConjTranspose
  };

  enum SideFlag {
    Right,
    Left
  };

  enum VectorFlag {
    Vectors,
    NoVectors
  };
}
