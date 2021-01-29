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
  pgesvd( const char* JOBU, const char* JOBVT, int64_t M, int64_t N,
         T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         T* S,
         T* U,  int64_t IU,  int64_t JU,  const scalapack_desc& DESCU, 
         T* VT, int64_t IVT, int64_t JVT, const scalapack_desc& DESCVT,
	 T* WORK, int64_t LWORK
  ); 




template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, int64_t>
  pgesvd( const char* JOBU, const char* JOBVT, int64_t M, int64_t N,
         T* A, int64_t IA, int64_t JA, const scalapack_desc& DESCA,
         detail::real_t<T>* S,
         T* U,  int64_t IU,  int64_t JU,  const scalapack_desc& DESCU, 
         T* VT, int64_t IVT, int64_t JVT, const scalapack_desc& DESCVT,
	 T* WORK, int64_t LWORK, detail::real_t<T>* RWORK
  ); 
           
}
}

