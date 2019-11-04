#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/descinit.hpp>

namespace scalapackpp {

BlockCyclicDist2D::BlockCyclicDist2D( const blacspp::Grid* grid,
  scalapack_int MB, scalapack_int NB, scalapack_int ISRC, scalapack_int JSRC
) : grid_(grid), mb_(MB), nb_(NB), isrc_(ISRC), jsrc_(JSRC){ }

BlockCyclicDist2D::BlockCyclicDist2D( const blacspp::Grid& grid,
  scalapack_int MB, scalapack_int NB, scalapack_int ISRC, scalapack_int JSRC
) : BlockCyclicDist2D( &grid, MB, NB, ISRC, JSRC ){ }

BlockCyclicDist2D::BlockCyclicDist2D( 
  const blacspp::Grid& grid, scalapack_int MB, scalapack_int NB
) : BlockCyclicDist2D( grid, MB, NB, 0, 0 ){ }

BlockCyclicDist2D::BlockCyclicDist2D() : BlockCyclicDist2D( nullptr, 0, 0, 0, 0 ){ }


BlockCyclicDist2D::BlockCyclicDist2D( const BlockCyclicDist2D& ) = default;
BlockCyclicDist2D::BlockCyclicDist2D( BlockCyclicDist2D&& ) noexcept = default;


bool BlockCyclicDist2D::is_valid() const {
  bool valid = grid_ and nb_ and mb_;
  if( valid ) 
    valid = valid and grid_->is_valid() and 
            (isrc_ < grid_->ipr()) and
            (jsrc_ < grid_->ipc());
  return valid;
}











std::pair< scalapack_int, scalapack_int > 
  BlockCyclicDist2D::get_local_dims( scalapack_int M, scalapack_int N ) {

  return scalapackpp::get_local_dims( *grid_, M, N, mb_, nb_, isrc_, jsrc_ );

}


std::pair< scalapack_desc, scalapack_int >
  BlockCyclicDist2D::descinit( scalapack_int M, scalapack_int N, scalapack_int LDD ){

  return scalapackpp::descinit( *grid_, M, N, mb_, nb_, isrc_, jsrc_, LDD );

}

scalapack_desc
  BlockCyclicDist2D::descinit_noerror( scalapack_int M, scalapack_int N, scalapack_int LDD ){

  return scalapackpp::descinit_noerror( *grid_, M, N, mb_, nb_, isrc_, jsrc_, LDD );

}



}
