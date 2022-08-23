#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "lace/vector.hpp"
#include <doctest/doctest.h>

consteval auto vector_op_nonconst() {
  lace::vector<double, 1> v{1.0};
  v[0] = 2.0;
  return v[0];
}

consteval auto vector_op_const() {
  const lace::vector<double, 1> v{1.0};
  return v[0];
}

TEST_CASE("lace/vector.hpp: vector::operator[]") {
  SUBCASE("nonconst") { CHECK(vector_op_nonconst() == 2.0); }
  SUBCASE("const") { CHECK(vector_op_const() == 1.0); }
}

TEST_CASE("lace/vector.hpp: vector::basis()") {
  auto v = lace::vector<double, 3>::basis(0);
  CHECK(v.at(0) == 1.0);
  CHECK(v.at(1) == 0.0);
  CHECK(v.at(2) == 0.0);
  v = lace::vector<double, 3>::basis(1);
  CHECK(v.at(0) == 0.0);
  CHECK(v.at(1) == 1.0);
  CHECK(v.at(2) == 0.0);
  v = lace::vector<double, 3>::basis(2);
  CHECK(v.at(0) == 0.0);
  CHECK(v.at(1) == 0.0);
  CHECK(v.at(2) == 1.0);
}

TEST_CASE("lace/vector.hpp: vector::to_array()") {
  const auto v = lace::vector<double, 3>{1.0, 2.0, 3.0}.to_array();
  CHECK(v[0] == 1.0);
  CHECK(v[1] == 2.0);
  CHECK(v[2] == 3.0);
}

TEST_CASE("lace/vector.hpp: vector::abs()") {
  const auto v = lace::vector<double, 2>{-1.0, 2.0}.abs();
  CHECK(v.at(0) == 1.0);
  CHECK(v.at(1) == 2.0);
}

TEST_CASE("lace/vector.hpp: vector::linorm()") {
  const auto n = lace::vector<double, 2>{-1.0, 2.0}.linorm();
  CHECK(n == 2.0);
}

TEST_CASE("lace/vector.hpp: vector::l1norm()") {
  const auto n = lace::vector<double, 2>{-1.0, 2.0}.l1norm();
  CHECK(n == 3.0);
}

TEST_CASE("lace/vector.hpp: vector::squared_l2norm()") {
  const auto n = lace::vector<double, 2>{-1.0, 2.0}.squared_l2norm();
  CHECK(n == 5.0);
}

TEST_CASE("lace/vector.hpp: vector::l2norm()") {
  const auto n =
      lace::vector<double, 2>{-1.0, 2.0}.l2norm() - std::sqrt(5.0) + 32.0;
  CHECK(n - 32.0 == 0.0);
}

TEST_CASE("lace/vector.hpp: vector::operator+(vector)") {
  auto v =
      lace::vector<double, 2>{-1.0, 2.0} + lace::vector<double, 2>{1.0, 1.0};
  CHECK(v.at(0) == 0.0);
  CHECK(v.at(1) == 3.0);
}

TEST_CASE("lace/vector.hpp: vector::operator+(double)") {
  auto v = lace::vector<double, 2>{-1.0, 2.0} + 1.0;
  CHECK(v.at(0) == 0.0);
  CHECK(v.at(1) == 3.0);
}

TEST_CASE("lace/vector.hpp: vector::operator-(vector)") {
  auto v =
      lace::vector<double, 2>{-1.0, 2.0} - lace::vector<double, 2>{1.0, 1.0};
  CHECK(v.at(0) == -2.0);
  CHECK(v.at(1) == 1.0);
}

TEST_CASE("lace/vector.hpp: vector::operator-(double)") {
  auto v = lace::vector<double, 2>{-1.0, 2.0} - 1.0;
  CHECK(v.at(0) == -2.0);
  CHECK(v.at(1) == 1.0);
}

TEST_CASE("lace/vector.hpp: vector::operator*(double)") {
  auto v = lace::vector<double, 2>{-1.0, 2.0} * 2.0;
  CHECK(v.at(0) == -2.0);
  CHECK(v.at(1) == 4.0);
}

TEST_CASE("lace/vector.hpp: vector::operator/(double)") {
  auto v = lace::vector<double, 2>{-1.0, 2.0} / 2.0;
  CHECK(v.at(0) == -0.5);
  CHECK(v.at(1) == 1.0);
}

TEST_CASE("lace/vector.hpp: vector::operator*(vector)") {
  auto d =
      lace::vector<double, 2>{-1.0, 2.0} * lace::vector<double, 2>{1.0, 1.0};
  CHECK(d == 1.0);
}
