/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/descinit.hpp>
#include <scalapackpp/util/type_conversions.hpp>

using scalapackpp::internal::scalapack_int;

extern "C" {
  void descinit_( scalapack_int* DESC, 
                  const scalapack_int* M, const scalapack_int* N, 
                  const scalapack_int* MB, const scalapack_int* NB,
                  const scalapack_int* ISRC, const scalapack_int* JSRC, 
                  const scalapack_int* ICONTEXT, const scalapack_int* LDD,
                  scalapack_int* INFO );

}

namespace scalapackpp::wrappers {

std::pair< scalapack_desc, int64_t > descinit( 
  int64_t M, int64_t N, int64_t MB, int64_t NB,
  int64_t ISRC, int64_t JSRC, int64_t ICONTEXT,
  int64_t LDD 
) {

  auto _M        = detail::to_scalapack_int( M    );
  auto _N        = detail::to_scalapack_int( N    );
  auto _MB       = detail::to_scalapack_int( MB   );
  auto _NB       = detail::to_scalapack_int( NB   );
  auto _ISRC     = detail::to_scalapack_int( ISRC );
  auto _JSRC     = detail::to_scalapack_int( JSRC );
  auto _ICONTEXT = detail::to_scalapack_int( ICONTEXT );

  scalapack_int LDD_use = std::max((int64_t)1, LDD);
  internal::scalapack_desc<scalapack_int> desc;

  scalapack_int INFO;
  descinit_( desc.data(), &_M, &_N, &_MB, &_NB, &_ISRC, &_JSRC, 
             &_ICONTEXT, &LDD_use, &INFO );

  return std::make_pair( detail::to_interface_int(desc), (int64_t)INFO );
}

}
