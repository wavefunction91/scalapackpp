#pragma once
#include <scalapackpp/types.hpp>
#include <string>

namespace scalapackpp::detail {

  extern std::string type_string( TransposeFlag );
  extern std::string type_string( SideFlag );
  extern std::string type_string( VectorFlag );

}
