/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/wrappers/descinit.hpp>

using scalapackpp::scalapack_int;

extern "C" {
  void descinit_( scalapack_int* DESC, 
                  const scalapack_int* M, const scalapack_int* N, 
                  const scalapack_int* MB, const scalapack_int* NB,
                  const scalapack_int* ISRC, const scalapack_int* JSRC, 
                  const scalapack_int* ICONTEXT, const scalapack_int* LDD,
                  scalapack_int* INFO );

}

namespace scalapackpp::wrappers {

std::pair< scalapack_desc, scalapack_int > descinit( 
  scalapack_int M, scalapack_int N, scalapack_int MB, scalapack_int NB,
  scalapack_int ISRC, scalapack_int JSRC, scalapack_int ICONTEXT,
  scalapack_int LDD 
) {


  scalapack_int LDD_use = std::max((scalapack_int)1, LDD);
  scalapack_desc desc;

  scalapack_int INFO;
  descinit_( desc.data(), &M, &N, &MB, &NB, &ISRC, &JSRC, &ICONTEXT, &LDD_use, 
             &INFO );

  return std::make_pair( desc, INFO );
}

}
