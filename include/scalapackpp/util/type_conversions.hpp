/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <blacspp/util/type_conversions.hpp>
#include <string>

namespace scalapackpp {
namespace detail      {

  extern std::string type_string( TransposeFlag );
  extern std::string type_string( SideFlag );
  extern std::string type_string( VectorFlag );
  extern std::string type_string( MatrixNorm );

  inline internal::scalapack_int to_scalapack_int( int64_t i ) {
    return blacspp::detail::to_blacs_int( i );
  }

  inline internal::scalapack_desc<internal::scalapack_int> 
    to_scalapack_int( const scalapack_desc& DESC ) {

    internal::scalapack_desc<internal::scalapack_int> _DESC;
    for( int64_t i = 0; i < internal::scalapack_desc_size; ++i )
      _DESC[i] = to_scalapack_int(DESC[i]);

    return _DESC;

  }

  inline scalapack_desc 
    to_interface_int( const internal::scalapack_desc<internal::scalapack_int> DESC ) {

    scalapack_desc _DESC;
    for( int64_t i = 0; i < internal::scalapack_desc_size; ++i )
      _DESC[i] = DESC[i];
   
    return _DESC;

  }

}
}
