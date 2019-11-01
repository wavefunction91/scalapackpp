#pragma once
#include <scalapackpp/wrappers/gemm.hpp>
#include <scalapackpp/util/type_conversions.hpp>

namespace scalapackpp {

template <typename T>
detail::enable_if_scalapack_supported_t<T>
  pgemm( TransposeFlag transa, TransposeFlag transb,
         scalapack_int M, scalapack_int N, scalapack_int K, T ALPHA, 
         const T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         const T* B, scalapack_int IB, scalapack_int JB, const scalapack_desc& DESCB,
         T BETA,
         T* C, scalapack_int IC, scalapack_int JC, const scalapack_desc& DESCC ) {

  auto TRANSA = detail::type_string( transa );
  auto TRANSB = detail::type_string( transb );
  
  wrappers::pgemm( TRANSA.c_str(), TRANSB.c_str(), M, N, K, ALPHA, A, IA, JA,
                   DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC ); 

}


}
