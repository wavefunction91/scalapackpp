/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/information.hpp>
#include <scalapackpp/descinit.hpp>

namespace scalapackpp {

BlockCyclicDist2D::BlockCyclicDist2D( const blacspp::Grid* grid,
  int64_t MB, int64_t NB, int64_t ISRC, int64_t JSRC
) : grid_(grid), mb_(MB), nb_(NB), isrc_(ISRC), jsrc_(JSRC){ }

BlockCyclicDist2D::BlockCyclicDist2D( const blacspp::Grid& grid,
  int64_t MB, int64_t NB, int64_t ISRC, int64_t JSRC
) : BlockCyclicDist2D( &grid, MB, NB, ISRC, JSRC ){ }

BlockCyclicDist2D::BlockCyclicDist2D( 
  const blacspp::Grid& grid, int64_t MB, int64_t NB
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











std::pair< int64_t, int64_t > 
  BlockCyclicDist2D::get_local_dims( int64_t M, int64_t N ) const {

  return scalapackpp::get_local_dims( *grid_, M, N, mb_, nb_, isrc_, jsrc_ );

}


std::pair< scalapack_desc, int64_t >
  BlockCyclicDist2D::descinit( int64_t M, int64_t N, int64_t LDD ) const {

  return scalapackpp::descinit( *grid_, M, N, mb_, nb_, isrc_, jsrc_, LDD );

}

scalapack_desc
  BlockCyclicDist2D::descinit_noerror( int64_t M, int64_t N, int64_t LDD ) const {

  return scalapackpp::descinit_noerror( *grid_, M, N, mb_, nb_, isrc_, jsrc_, LDD );

}



}
