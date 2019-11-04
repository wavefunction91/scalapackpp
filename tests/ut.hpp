#include <catch2/catch.hpp>
#include <random>
#include <vector>
#include <iostream>

#define SCALAPACKPP_REAL_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, float, double )

#define SCALAPACKPP_COMPLEX_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, scalapackpp::scomplex, scalapackpp::dcomplex)

#define SCALAPACKPP_TEST_CASE(NAME, CAT)\
TEMPLATE_TEST_CASE(NAME,CAT, float, double, scalapackpp::scomplex, scalapackpp::dcomplex)

