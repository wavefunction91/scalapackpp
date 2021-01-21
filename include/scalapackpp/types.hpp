/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once

#include <blacspp/types.hpp>
#include <array>

namespace scalapackpp {

namespace internal {
  using scalapack_int = blacspp::internal::blacs_int;
  using scomplex      = blacspp::internal::scomplex;
  using dcomplex      = blacspp::internal::dcomplex;

  static constexpr size_t scalapack_desc_size = 9;

  static constexpr size_t _DTYPE_A = 0;
  static constexpr size_t _CTXT_A  = 1;
  static constexpr size_t _M_A     = 2;
  static constexpr size_t _N_A     = 3;
  static constexpr size_t _MB_A    = 4;
  static constexpr size_t _NB_A    = 5;
  static constexpr size_t _RSRC_A  = 6;
  static constexpr size_t _CSRC_A  = 7;
  static constexpr size_t _LLD_A   = 8;


  template <typename Integral>
  using scalapack_desc = std::array< Integral, scalapack_desc_size >;
}

  using scalapack_desc = internal::scalapack_desc<int64_t>;


  struct block_cyclic_coordinate {
    int64_t process_id;
    int64_t block_idx;
    int64_t local_idx;
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


  enum MatrixNorm {
    FrobeniusNorm,
    InfinityNorm,
    OneNorm,
    AbsMax
  };
}
