/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/util/type_traits.hpp>

namespace scalapackpp {
namespace wrappers    {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  psyev( const char* JOBZ, const char* UPLO, int64_t N,
         T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         T* W,
         T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ,
         T* WORK, int64_t LWORK );

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, int64_t>
  psyevd( const char* JOBZ, const char* UPLO, int64_t N,
         T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         T* W,
         T* Z, int64_t IZ, int64_t JZ, const scalapack_desc& DESCZ,
         T* WORK, int64_t LWORK, internal::scalapack_int* IWORK, int64_t LIWORK );

}
}
