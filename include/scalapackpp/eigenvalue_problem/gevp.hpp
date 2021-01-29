/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/pblas/trsm.hpp>
#include <scalapackpp/factorizations/potrf.hpp>
#include <scalapackpp/eigenvalue_problem/sevp.hpp>

#include <scalapackpp/wrappers/eigenvalue_problem/sygst.hpp>
#include <scalapackpp/wrappers/eigenvalue_problem/hegst.hpp>

namespace scalapackpp {

template <
  typename T,
  detail::enable_if_scalapack_real_supported_t<T,bool> = true
>
std::pair<int64_t, T>
  gen_to_std_evp( int64_t IBTYPE, blacspp::Triangle uplo, int64_t N,
                        T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
                  const T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {

  auto UPLO = blacspp::detail::type_string(uplo);
  T SCALE;
  auto info = 
    wrappers::psygst( IBTYPE, UPLO.c_str(), N, A, IA, JA, DESCA, B, IB, JB, DESCB, &SCALE );

  return std::make_pair( info, SCALE );

}

template <
  typename T,
  detail::enable_if_scalapack_complex_supported_t<T,bool> = true
>
std::pair<int64_t, detail::real_t<T>>
  gen_to_std_evp( int64_t IBTYPE, blacspp::Triangle uplo, int64_t N,
                        T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
                  const T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB ) {

  auto UPLO = blacspp::detail::type_string(uplo);
  detail::real_t<T> SCALE;
  auto info = 
    wrappers::phegst( IBTYPE, UPLO.c_str(), N, A, IA, JA, DESCA, B, IB, JB, DESCB, &SCALE );

  return std::make_pair( info, SCALE );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  hereig_gen( VectorFlag jobz, blacspp::Triangle uplo, int64_t N,
              T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
              T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB,
              detail::real_t<T>* W,
              T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ )
  {

    // Perform a Cholesky factorization on B
    auto info = ppotrf( uplo, N, B, IB, JB, DESCB );
    if( info ) return info;

#if 0
    // Transform A -> L**-1 * A * L**-H
    ptrsm( SideFlag::Left, uplo, TransposeFlag::NoTranspose, 
           blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
           A, IA, JA, DESCA );
    ptrsm( SideFlag::Right, uplo, TransposeFlag::ConjTranspose, 
           blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
           A, IA, JA, DESCA );
#else
    detail::real_t<T> SCALE;
    std::tie( info, SCALE) = 
      gen_to_std_evp( 1, uplo, N, A, IA, JA, DESCA, B, IB, JB, DESCB );

    // Check that scale is 1
    assert( std::abs( SCALE - 1 ) < std::numeric_limits<detail::real_t<T>>::epsilon() );
#endif

    // Solve SEVP
    info = hereig( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );
    if( info ) return info;


    // If eigenvectors requested, backtransform Z
    if( jobz == VectorFlag::Vectors ) {
      if( uplo == blacspp::Triangle::Lower )
        ptrsm( SideFlag::Left, uplo, TransposeFlag::ConjTranspose,
               blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
               Z, IZ, JZ, DESCZ );
      else
        ptrsm( SideFlag::Left, uplo, TransposeFlag::NoTranspose,
               blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
               Z, IZ, JZ, DESCZ );
    }
     
    return info;
}
             
template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  hereigd_gen( VectorFlag jobz, blacspp::Triangle uplo, int64_t N,
              T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
              T* B, int64_t IB, int64_t JB, const scalapack_desc& DESCB,
              detail::real_t<T>* W,
              T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ )
  {

    // Perform a Cholesky factorization on B
    auto info = ppotrf( uplo, N, B, IB, JB, DESCB );
    if( info ) return info;

#if 0
    // Transform A -> L**-1 * A * L**-H
    ptrsm( SideFlag::Left, uplo, TransposeFlag::NoTranspose, 
           blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
           A, IA, JA, DESCA );
    ptrsm( SideFlag::Right, uplo, TransposeFlag::ConjTranspose, 
           blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
           A, IA, JA, DESCA );
#else
    detail::real_t<T> SCALE;
    std::tie( info, SCALE) = 
      gen_to_std_evp( 1, uplo, N, A, IA, JA, DESCA, B, IB, JB, DESCB );

    // Check that scale is 1
    assert( std::abs( SCALE - 1 ) < std::numeric_limits<detail::real_t<T>>::epsilon() );
#endif

    // Solve SEVP
    info = hereigd( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );
    if( info ) return info;


    // If eigenvectors requested, backtransform Z
    if( jobz == VectorFlag::Vectors ) {
      if( uplo == blacspp::Triangle::Lower )
        ptrsm( SideFlag::Left, uplo, TransposeFlag::ConjTranspose,
               blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
               Z, IZ, JZ, DESCZ );
      else
        ptrsm( SideFlag::Left, uplo, TransposeFlag::NoTranspose,
               blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
               Z, IZ, JZ, DESCZ );
    }
     
    return info;
}

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  hereig_gen( VectorFlag jobz, blacspp::Triangle uplo, 
              BlockCyclicMatrix<T>& A, 
              BlockCyclicMatrix<T>& B, 
              detail::real_t<T>* W, BlockCyclicMatrix<T>& Z ) {

  // TODO Sanity check
  return hereig_gen( jobz, uplo, A.m(), A.data(), 1, 1, A.desc(), 
                     B.data(), 1, 1, B.desc(), W, Z.data(), 1, 1,
                     Z.desc() );

}

template <typename T>
detail::enable_if_scalapack_supported_t<T, int64_t>
  hereigd_gen( VectorFlag jobz, blacspp::Triangle uplo, 
               BlockCyclicMatrix<T>& A, 
               BlockCyclicMatrix<T>& B, 
               detail::real_t<T>* W, BlockCyclicMatrix<T>& Z ) {

  // TODO Sanity check
  return hereigd_gen( jobz, uplo, A.m(), A.data(), 1, 1, A.desc(), 
                      B.data(), 1, 1, B.desc(), W, Z.data(), 1, 1,
                      Z.desc() );

}


}
