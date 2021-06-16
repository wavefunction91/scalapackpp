/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/util/type_traits.hpp>
#include <vector>

namespace scalapackpp {

template <typename T>
class BlockCyclicMatrix {

  using Grid = blacspp::Grid;

  int64_t M_;
  int64_t N_;

  int64_t M_local_;
  int64_t N_local_;

  BlockCyclicDist2D dist_;
  std::vector<T>    local_data_;
  scalapack_desc    desc_;

public:

  constexpr BlockCyclicMatrix() noexcept : M_(0), N_(0) { };

  BlockCyclicMatrix( const Grid& grid, int64_t M, int64_t N, int64_t MB,
                     int64_t NB, int64_t ISRC, int64_t JSRC ) :
    M_(M), N_(N), dist_( grid, MB, NB, ISRC, JSRC ) {

    std::tie( M_local_, N_local_ ) = dist_.get_local_dims( M, N );
    
    int64_t desc_info;
    std::tie( desc_, desc_info ) = dist_.descinit( M, N, M_local_ );
    if( desc_info )
      throw std::runtime_error("DESCINIT FAILED");

    local_data_.resize( M_local_ * N_local_ );

  }


  BlockCyclicMatrix( const Grid& grid, int64_t M, int64_t N, int64_t MB,
                     int64_t NB ) :
    BlockCyclicMatrix( grid, M, N, MB, NB, 0, 0 ) { }


  BlockCyclicMatrix( const BlockCyclicMatrix& )     = default;
  BlockCyclicMatrix( BlockCyclicMatrix&& ) noexcept = default;

  BlockCyclicMatrix& operator=( const BlockCyclicMatrix& )     = default;
  BlockCyclicMatrix& operator=( BlockCyclicMatrix&& ) noexcept = default;



  inline auto m()       const { return M_;          }
  inline auto n()       const { return N_;          }
  inline auto m_local() const { return M_local_;    }
  inline auto n_local() const { return N_local_;    }
  inline auto mb()      const { return dist_.mb();  }
  inline auto nb()      const { return dist_.nb();  }
  inline auto isrc()    const { return dist_.isrc(); }
  inline auto jsrc()    const { return dist_.jsrc(); }

  inline auto& desc() const { return desc_; }
  inline auto  data() const { return local_data_.data(); }
  inline auto  data()       { return local_data_.data(); }
  inline auto  local_size() const { return local_data_.size(); }

  inline auto& dist() const { return dist_; }


  inline void scatter_to( int64_t M, int64_t N, const T* A, int64_t LDA, 
                          int64_t ISRC, int64_t JSRC ) {
    // TODO: Check M/N
    dist_.scatter( M_, N_, A, LDA, data(), M_local_, ISRC, JSRC );

  }

  inline void gather_from( int64_t M, int64_t N, T* A, int64_t LDA, 
                           int64_t IDEST, int64_t JDEST ) const {
    // TODO: Check M/N
    dist_.gather( M_, N_, A, LDA, data(), M_local_, IDEST, JDEST );

  }



  inline auto begin()        { return local_data_.begin(); };
  inline auto end()          { return local_data_.end();   };
  inline auto begin()  const { return local_data_.begin(); };
  inline auto end()    const { return local_data_.end();   };
  inline auto cbegin() const { return local_data_.cbegin(); };
  inline auto cend()   const { return local_data_.cend();   };

};



template <typename T>
inline detail::enable_if_scalapack_supported_t<T> 
  redistribute( const BlockCyclicMatrix<T>& A, BlockCyclicMatrix<T>& B ) {

  // TODO Sanity check on A/B
  // TODO Sanity check / handle union of contexts
  wrappers::pgemr2d( A.m(), A.n(), A.data(), 1, 1, A.desc(), B.data(), 1, 1, 
                     B.desc(), A.context() );

}

}
