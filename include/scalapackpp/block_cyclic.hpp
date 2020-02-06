#pragma once
#include <scalapackpp/types.hpp>
#include <blacspp/grid.hpp>
#include <scalapackpp/scatter_gather.hpp>

namespace scalapackpp {

class BlockCyclicDist2D {

  using Grid = blacspp::Grid;

  const Grid*   grid_ = nullptr;
  scalapack_int mb_;
  scalapack_int nb_;
  scalapack_int isrc_;
  scalapack_int jsrc_;


  BlockCyclicDist2D( const Grid* grid, scalapack_int MB, scalapack_int NB,
                     scalapack_int ISRC, scalapack_int JSRC);

public:

  bool is_valid() const;

  BlockCyclicDist2D();

  BlockCyclicDist2D( const Grid& grid, scalapack_int MB, scalapack_int NB,
                     scalapack_int ISRC, scalapack_int JSRC);

  BlockCyclicDist2D( const Grid& grid, scalapack_int MB, scalapack_int NB );


  BlockCyclicDist2D( const BlockCyclicDist2D& );
  BlockCyclicDist2D( BlockCyclicDist2D&& ) noexcept ;

  inline auto mb() const { return mb_; }
  inline auto nb() const { return nb_; }

  std::pair< scalapack_int, scalapack_int > 
    get_local_dims( scalapack_int M, scalapack_int N );

  std::pair< scalapack_desc, scalapack_int >
    descinit( scalapack_int M, scalapack_int N, scalapack_int LDD );

  scalapack_desc
    descinit_noerror( scalapack_int M, scalapack_int N, scalapack_int LDD );


  template <typename T>
  void scatter( scalapack_int M, scalapack_int N, const T* A, scalapack_int LDA,
                T* A_local, scalapack_int LDA_local,
                scalapack_int ISRC, scalapack_int JSRC ) const {

    scalapackpp::scatter( *grid_, M, N, mb_, nb_, A, LDA, ISRC, JSRC,
                          A_local, LDA_local, isrc_, jsrc_ );

  }


  template <typename T>
  void gather( scalapack_int M, scalapack_int N, T* A, scalapack_int LDA,
               const T* A_local, scalapack_int LDA_local,
               scalapack_int IDEST, scalapack_int JDEST ) const {

    scalapackpp::gather( *grid_, M, N, mb_, nb_, A, LDA, IDEST, JDEST,
                          A_local, LDA_local, isrc_, jsrc_ );

  }



  inline std::pair< scalapack_int, scalapack_int >
    owner_coordinate( scalapack_int I, scalapack_int J ) const noexcept {

    return { (I / mb_) % grid_->npr(), (J / nb_) % grid_->npc() };

  }



  inline bool i_own( scalapack_int I, scalapack_int J ) const noexcept {
    auto [pr, pc] = owner_coordinate( I, J );
    return grid_->ipr() == pr and grid_->ipc() == pc;
  }

  inline std::pair< scalapack_int, scalapack_int >
    local_indx( scalapack_int I, scalapack_int J ) const noexcept {

    auto l = I / (mb_ * grid_->npr());
    auto m = J / (nb_ * grid_->npc());

    return { l * mb_ + (I % mb_), m * nb_ + (J % nb_) };

  }


  inline const blacspp::Grid& grid() const { return *grid_; }

};

}
