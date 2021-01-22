/**
 *  This file is a part of scalapackpp (see LICENSE)
 *
 *  Copyright (c) 2019-2020 David Williams-Young
 *  All rights reserved
 */
#pragma once
#include <scalapackpp/types.hpp>
#include <type_traits>
#include <cassert>

#if __cplusplus < 201402L
namespace std {
template <typename T, typename U>
using enable_if_t = typename enable_if<T,U>::type;
}
#endif

namespace scalapackpp {
namespace detail {

template <typename T>
struct scalapack_real_supported : public std::false_type { };
template <typename T>
struct scalapack_complex_supported : public std::false_type { };

template <>
struct scalapack_real_supported< float >    : public std::true_type { };
template <>
struct scalapack_real_supported< double >   : public std::true_type { };

template <>
struct scalapack_complex_supported< internal::scomplex > : public std::true_type { };
template <>
struct scalapack_complex_supported< internal::dcomplex > : public std::true_type { };


template <typename T>
struct scalapack_supported {
  static constexpr bool value = scalapack_real_supported<T>::value or scalapack_complex_supported<T>::value;
};


#if __cplusplus >= 201703L

template <typename T>
inline constexpr bool scalapack_real_supported_v = 
  scalapack_real_supported<T>::value;

template <typename T>
inline constexpr bool scalapack_complex_supported_v = 
  scalapack_complex_supported<T>::value;

template <typename T>
inline constexpr bool scalapack_supported_v = 
  scalapack_supported<T>::value;

#endif


template <typename T, typename U = void>
using enable_if_scalapack_real_supported_t = 
  typename std::enable_if< scalapack_real_supported<T>::value, U >::type;

template <typename T, typename U = void>
using enable_if_scalapack_complex_supported_t = 
  typename std::enable_if< scalapack_complex_supported<T>::value, U >::type;

template <typename T, typename U = void>
using enable_if_scalapack_supported_t = 
  typename std::enable_if< scalapack_supported<T>::value, U >::type;




template <typename T>
struct real {
  using type = T;
};

template <typename T>
struct real< std::complex<T> > {
  using type = T;
};

template <typename T>
using real_t = typename real<T>::type;



}
}
