#pragma once
#include <scalapackpp/types.hpp>
#include <blacspp/grid.hpp>

namespace scalapackpp {

  scalapack_int numroc( scalapack_int N, scalapack_int NB,
    scalapack_int IPROC, scalapack_int ISRC, 
    scalapack_int NPROC );

  std::pair<scalapack_int, scalapack_int>
    get_local_dims( const blacspp::Grid& grid,
      scalapack_int M, scalapack_int N,
      scalapack_int MB, scalapack_int NB,
      scalapack_int ISRC, scalapack_int JSRC );

/*
  std::pair< block_cyclic_coordinate, block_cyclic_coordinate >
    local_index( const blacspp::Grid& grid,
      scalapack_int MB, scalapack_int NB,
      scalapack_int I, scalapack_int J
    );

  std::pair< scalapack_int, scalapack_int >
    global_index( const blacspp::Grid& grid,
      scalapack_int MB, scalapack_int NB,
      scalapack_int I, scalapack_int J
    );
*/
}
