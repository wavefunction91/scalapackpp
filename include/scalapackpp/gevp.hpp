#pragma once
#include <scalapackpp/sevp.hpp>
#include <scalapackpp/trsm.hpp>
#include <scalapackpp/potrf.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  hereig_gen( VectorFlag jobz, blacspp::Triangle uplo, scalapack_int N,
              T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
              T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
              detail::real_t<T>* W,
              T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ )
  {

    // Perform a Cholesky factorization on B
    auto info = ppotrf( uplo, N, B, IB, JB, DESCB );
    if( info ) return info;

    // Transform A -> L**-1 * A * L**-H
    ptrsm( SideFlag::Left, uplo, TransposeFlag::NoTranspose, 
           blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
           A, IA, JA, DESCA );
    ptrsm( SideFlag::Right, uplo, TransposeFlag::ConjTranspose, 
           blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
           A, IA, JA, DESCA );

    // Solve SEVP
    info = hereig( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );
    if( info ) return info;


    // If eigenvectors requested, backtransform Z
    if( jobz == VectorFlag::Vectors )
      ptrsm( SideFlag::Left, uplo, TransposeFlag::ConjTranspose,
             blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
             Z, IZ, JZ, DESCZ );
     
    return info;
}
             
template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  hereigd_gen( VectorFlag jobz, blacspp::Triangle uplo, scalapack_int N,
              T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
              T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
              detail::real_t<T>* W,
              T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ )
  {

    // Perform a Cholesky factorization on B
    auto info = ppotrf( uplo, N, B, IB, JB, DESCB );
    if( info ) return info;

    // Transform A -> L**-1 * A * L**-H
    ptrsm( SideFlag::Left, uplo, TransposeFlag::NoTranspose, 
           blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
           A, IA, JA, DESCA );
    ptrsm( SideFlag::Right, uplo, TransposeFlag::ConjTranspose, 
           blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
           A, IA, JA, DESCA );

    // Solve SEVP
    info = hereigd( jobz, uplo, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ );
    if( info ) return info;


    // If eigenvectors requested, backtransform Z
    if( jobz == VectorFlag::Vectors )
      ptrsm( SideFlag::Left, uplo, TransposeFlag::ConjTranspose,
             blacspp::Diagonal::NonUnit, N, N, 1., B, IB, JB, DESCB,
             Z, IZ, JZ, DESCZ );
     
    return info;
}


}
