#pragma once
#include <scalapackpp/types.hpp>

namespace scalapackpp::wrappers {

std::pair< scalapack_desc, scalapack_int > descinit( 
  scalapack_int M, scalapack_int N, scalapack_int MB, scalapack_int NB,
  scalapack_int ISRC, scalapack_int JSRC, scalapack_int ICONTEXT,
  scalapack_int LDD 
);

}
