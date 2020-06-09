/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <scalapackpp/util/sfinae.hpp>

namespace scalapackpp::wrappers {

template <typename T>
detail::enable_if_scalapack_real_supported_t<T, scalapack_int>
  pgesvd( const char* JOBU, const char* JOBVT, scalapack_int M, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         T* S,
         T* U,  scalapack_int IU,  scalapack_int JU,  const scalapack_desc& DESCU, 
         T* VT, scalapack_int IVT, scalapack_int JVT, const scalapack_desc& DESCVT,
	 T* WORK, scalapack_int LWORK
  ); 




template <typename T>
detail::enable_if_scalapack_complex_supported_t<T, scalapack_int>
  pgesvd( const char* JOBU, const char* JOBVT, scalapack_int M, scalapack_int N,
         T* A, scalapack_int IA, scalapack_int JA, const scalapack_desc& DESCA,
         detail::real_t<T>* S,
         T* U,  scalapack_int IU,  scalapack_int JU,  const scalapack_desc& DESCU, 
         T* VT, scalapack_int IVT, scalapack_int JVT, const scalapack_desc& DESCVT,
	 T* WORK, scalapack_int LWORK, detail::real_t<T>* RWORK
  ); 
           
}

