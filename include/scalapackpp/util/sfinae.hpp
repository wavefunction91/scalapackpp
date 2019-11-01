#pragma once
#include <scalapackpp/types.hpp>
#include <type_traits>

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


template <typename T, typename U = void>
using enable_if_scalapack_real_supported_t = 
  typename std::enable_if< scalapack_real_supported<T>::value, U >::type;

template <typename T, typename U = void>
using enable_if_scalapack_complex_supported_t = 
  typename std::enable_if< scalapack_complex_supported<T>::value, U >::type;

template <typename T, typename U = void>
using enable_if_scalapack_supported_t = 
  typename std::enable_if< 
    scalapack_real_supported<T>::value or scalapack_complex_supported<T>::value, U 
  >::type;




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
