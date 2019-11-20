#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/util/sfinae.hpp>

namespace scalapackpp::wrappers {

template <typename T>
detail::enable_if_scalapack_supported_t<T, scalapack_int>
  pgetri( scalapack_int N, T* A, scalapack_int IA, scalapack_int JA, 
          const scalapack_desc& DESCA, const scalapack_int* IPIV, 
          T* WORK, scalapack_int LWORK, scalapack_int* IWORK, scalapack_int LIWORK );


}




