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


namespace scalapackpp {
namespace detail {

// Define local enable_if_t for C++11 compilation
#if __cplusplus < 201402L
  // < C++14
  template <bool B, typename T = void>
  using enable_if_t = typename std::enable_if<B,T>::type;
#else
  // >= C++14
  template <bool B, typename T = void>
  using enable_if_t = std::enable_if_t<B,T>;
#endif

// Define local type_dentity for C++11/14/17
#if __cplusplus > 201703L
  // >= C++20
  template <typename T>
  using type_identity = std::type_identity<T>;

  template <typename T>
  using type_identity_t = std::type_identity_t<T>;
#else
  // < C++20
  template <typename T>
  struct type_identity {
    using type = T;
  };

  template <typename T>
  using type_identity_t = typename type_identity<T>::type;
#endif


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


// Define local _v's for >= C++17
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
