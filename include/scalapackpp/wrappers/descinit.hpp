/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>

namespace scalapackpp {
namespace wrappers    {

std::pair< scalapack_desc, int64_t > descinit( 
  int64_t M, int64_t N, int64_t MB, int64_t NB,
  int64_t ISRC, int64_t JSRC, int64_t ICONTEXT,
  int64_t LDD 
);

}
}
