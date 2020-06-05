/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/descinit.hpp>
#include <scalapackpp/wrappers/gemr2d.hpp>
#include <blacspp/grid.hpp>
#include <tuple>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  scatter( const blacspp::Grid& grid, 
    scalapack_int M, scalapack_int N, scalapack_int MB, scalapack_int NB,
    const T* A, scalapack_int LDA, scalapack_int ILOC, scalapack_int JLOC,  
    T* A_local, scalapack_int LDLOCA, scalapack_int ISRC, scalapack_int JSRC ) {

  auto desc_a     = descinit_noerror( grid, M, N, M,  N,  ILOC, JLOC, LDA    );
  auto desc_loc_a = descinit_noerror( grid, M, N, MB, NB, ISRC, JSRC, LDLOCA );

  auto context = grid.context();
  wrappers::pgemr2d( M, N, A, 1, 1, desc_a, A_local, 1, 1, desc_loc_a, context );

}
    
template <typename T>
detail::enable_if_scalapack_supported_t<T>
  gather( const blacspp::Grid& grid, 
    scalapack_int M, scalapack_int N, scalapack_int MB, scalapack_int NB,
    T* A, scalapack_int LDA, scalapack_int ILOC, scalapack_int JLOC,  
    const T* A_local, scalapack_int LDLOCA, scalapack_int ISRC, scalapack_int JSRC 
  ) {

  auto desc_a     = descinit_noerror( grid, M, N, M,  N,  ILOC, JLOC, LDA    );
  auto desc_loc_a = descinit_noerror( grid, M, N, MB, NB, ISRC, JSRC, LDLOCA );

  auto context = grid.context();
  wrappers::pgemr2d( M, N, A_local, 1, 1, desc_loc_a, A, 1, 1, desc_a, context );

}
    

}
