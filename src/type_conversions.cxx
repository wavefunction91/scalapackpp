#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp::detail {

std::string type_string( TransposeFlag trans ) {
  if( trans == NoTranspose )     return std::string( "N" );
  else if( trans == Transpose )  return std::string( "T" );
  else                           return std::string( "C" );
}

}
