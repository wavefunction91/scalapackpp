#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/util/sfinae.hpp>

namespace scalapackpp::wrappers {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, scalapack_int>
  psyev( const char* JOBZ, const char* UPLO, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         T* W,
         T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
         T* WORK, scalapack_int LWORK );

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, scalapack_int>
  psyevd( const char* JOBZ, const char* UPLO, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         T* W,
         T* Z, scalapack_int IZ, scalapack_int JZ, const scalapack_desc& DESCZ,
         T* WORK, scalapack_int LWORK, scalapack_int* IWORK, scalapack_int LIWORK );

}
