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

namespace scalapackpp {

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
