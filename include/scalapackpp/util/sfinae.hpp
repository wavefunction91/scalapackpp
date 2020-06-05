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

namespace scalapackpp::detail {

template <typename T>
struct scalapack_real_supported : public std::false_type { };
template <typename T>
struct scalapack_complex_supported : public std::false_type { };

template <>
struct scalapack_real_supported< float >    : public std::true_type { };
template <>
struct scalapack_real_supported< double >   : public std::true_type { };

template <>
struct scalapack_complex_supported< scomplex > : public std::true_type { };
template <>
struct scalapack_complex_supported< dcomplex > : public std::true_type { };

template <typename T>
inline constexpr bool scalapack_real_supported_v = 
  scalapack_real_supported<T>::value;

template <typename T>
inline constexpr bool scalapack_complex_supported_v = 
  scalapack_complex_supported<T>::value;

template <typename T>
inline constexpr bool scalapack_supported_v = 
  scalapack_real_supported_v<T> or scalapack_complex_supported_v<T>;


template <typename T, typename U = void>
using enable_if_scalapack_real_supported_t = 
  typename std::enable_if< scalapack_real_supported_v<T>, U >::type;

template <typename T, typename U = void>
using enable_if_scalapack_complex_supported_t = 
  typename std::enable_if< scalapack_complex_supported_v<T>, U >::type;

template <typename T, typename U = void>
using enable_if_scalapack_supported_t = 
  typename std::enable_if< scalapack_supported_v<T>, U >::type;




template <typename T>
struct real {
  using type = T;
};

template<>
struct real< scomplex > {
  using type = float;
};
template<>
struct real< dcomplex > {
  using type = double;
};

template <typename T>
using real_t = typename real<T>::type;



}
