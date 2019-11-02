#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp::detail {

std::string type_string( TransposeFlag trans ) {
  if( trans == NoTranspose )     return std::string( "N" );
  else if( trans == Transpose )  return std::string( "T" );
  else                           return std::string( "C" );
}

std::string type_string( SideFlag side ) {
  if( side == Right ) return std::string( "R" );
  else                return std::string( "L" );
}

}
