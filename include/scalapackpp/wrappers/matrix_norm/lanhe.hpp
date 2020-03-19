#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/util/sfinae.hpp>

namespace scalapackpp::wrappers {

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, detail::real_t<T>>
  planhe( const char* NORM, const char* UPLO, scalapack_int N, const T* A, 
          scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
          detail::real_t<T>* WORK );

}

