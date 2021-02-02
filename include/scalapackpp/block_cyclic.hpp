/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <blacspp/grid.hpp>
#include <scalapackpp/scatter_gather.hpp>

#include <tuple>

namespace scalapackpp {

class BlockCyclicDist2D {

  using Grid = blacspp::Grid;

  Grid grid_;
  int64_t mb_;
  int64_t nb_;
  int64_t isrc_;
  int64_t jsrc_;

public:

  bool is_valid() const;

  BlockCyclicDist2D();

  BlockCyclicDist2D( const Grid& grid, int64_t MB, int64_t NB,
                     int64_t ISRC, int64_t JSRC);

  BlockCyclicDist2D( const Grid& grid, int64_t MB, int64_t NB );


  BlockCyclicDist2D( const BlockCyclicDist2D& );
  BlockCyclicDist2D( BlockCyclicDist2D&& ) noexcept ;

  inline auto mb() const { return mb_; }
  inline auto nb() const { return nb_; }
  inline auto isrc() const { return isrc_; }
  inline auto jsrc() const { return jsrc_; }
  

  std::pair< int64_t, int64_t > 
    get_local_dims( int64_t M, int64_t N ) const;

  std::pair< scalapack_desc, int64_t >
    descinit( int64_t M, int64_t N, int64_t LDD ) const;

  scalapack_desc
    descinit_noerror( int64_t M, int64_t N, int64_t LDD ) const;


  template <typename T>
  void scatter( int64_t M, int64_t N, const T* A, int64_t LDA,
                T* A_local, int64_t LDA_local,
                int64_t ISRC, int64_t JSRC ) const {

    scalapackpp::scatter( grid_, M, N, mb_, nb_, A, LDA, ISRC, JSRC,
                          A_local, LDA_local, isrc_, jsrc_ );

  }


  template <typename T>
  void gather( int64_t M, int64_t N, T* A, int64_t LDA,
               const T* A_local, int64_t LDA_local,
               int64_t IDEST, int64_t JDEST ) const {

    scalapackpp::gather( grid_, M, N, mb_, nb_, A, LDA, IDEST, JDEST,
                          A_local, LDA_local, isrc_, jsrc_ );

  }



  inline std::pair< int64_t, int64_t >
    owner_coordinate( int64_t I, int64_t J ) const noexcept {

    return { (I / mb_) % grid_.npr(), (J / nb_) % grid_.npc() };

  }



  inline bool i_own( int64_t I, int64_t J ) const noexcept {
    int64_t pr,pc;
    std::tie( pr, pc ) = owner_coordinate( I, J );
    return grid_.ipr() == pr and grid_.ipc() == pc;
  }

  inline std::pair< int64_t, int64_t >
    local_indx( int64_t I, int64_t J ) const noexcept {

    auto l = I / (mb_ * grid_.npr());
    auto m = J / (nb_ * grid_.npc());

    return { l * mb_ + (I % mb_), m * nb_ + (J % nb_) };

  }


  inline const blacspp::Grid& grid() const { return grid_; }

};

}
