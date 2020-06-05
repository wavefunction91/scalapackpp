/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <string>

namespace scalapackpp::detail {

  extern std::string type_string( TransposeFlag );
  extern std::string type_string( SideFlag );
  extern std::string type_string( VectorFlag );
  extern std::string type_string( MatrixNorm );

}
