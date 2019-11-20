#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/util/sfinae.hpp>

namespace scalapackpp::wrappers {

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, scalapack_int>
  pheev( const char* JOBZ, const char* UPLO, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         detail::real_t<T>* W,
         T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
         T* WORK, scalapack_int LWORK, 
         detail::real_t<T>* RWORK, scalapack_int LRWORK );

template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, scalapack_int>
  pheevd( const char* JOBZ, const char* UPLO, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         detail::real_t<T>* W,
         T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
         T* WORK, scalapack_int LWORK, 
         detail::real_t<T>* RWORK, scalapack_int LRWORK, 
         scalapack_int* IWORK, scalapack_int LIWORK );

}
