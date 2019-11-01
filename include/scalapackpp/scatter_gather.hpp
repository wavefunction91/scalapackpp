#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/wrappers/gemr2d.hpp>
#include <scalapackpp/wrappers/descinit.hpp>
#include <blacspp/grid.hpp>
#include <tuple>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  scatter( const blacspp::Grid& grid, 
    scalapack_int M, scalapack_int N, scalapack_int MB, scalapack_int NB,
    const T* A, scalapack_int LDA, scalapack_int ILOC, scalapack_int JLOC,  
    T* A_local, scalapack_int LDLOCA, scalapack_int ISRC, scalapack_int JSRC ) {

  auto context = grid.context();
  auto [desc_a    ,i1] = wrappers::descinit( M, N, M,  N,  ILOC, JLOC, context, LDA    );
  auto [desc_loc_a,i2] = wrappers::descinit( M, N, MB, NB, ISRC, JSRC, context, LDLOCA );

  wrappers::pgemr2d( M, N, A, 1, 1, desc_a, A_local, 1, 1, desc_loc_a, context );

}
    
template <typename T>
detail::enable_if_scalapack_supported_t<T>
  gather( const blacspp::Grid& grid, 
    scalapack_int M, scalapack_int N, scalapack_int MB, scalapack_int NB,
    T* A, scalapack_int LDA, scalapack_int ILOC, scalapack_int JLOC,  
    const T* A_local, scalapack_int LDLOCA, scalapack_int ISRC, scalapack_int JSRC 
  ) {

  auto context = grid.context();
  auto [desc_a    ,i1] = wrappers::descinit( M, N, M,  N,  ILOC, JLOC, context, LDA    );
  auto [desc_loc_a,i2] = wrappers::descinit( M, N, MB, NB, ISRC, JSRC, context, LDLOCA );

  wrappers::pgemr2d( M, N, A_local, 1, 1, desc_loc_a, A, 1, 1, desc_a, context );

}
    

}
