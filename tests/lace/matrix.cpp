#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "lace/matrix.hpp"
#include <doctest/doctest.h>

consteval auto matrix_op_nonconst() {
  lace::square_matrix<double, 1> m{1.0};
  m[0][0] = 2.0;
  return m[0][0];
}

consteval auto matrix_op_const() {
  const lace::square_matrix<double, 1> m{1.0};
  return m[0][0];
}

TEST_CASE("lace/matrix.hpp: square_matrix::operator[]") {
  SUBCASE("nonconst") { CHECK(matrix_op_nonconst() == 2.0); }
  SUBCASE("const") { CHECK(matrix_op_const() == 1.0); }
}

TEST_CASE("lace/matrix.hpp: square_matrix::identity()") {
  auto m = lace::square_matrix<double, 2>::identity();
  CHECK(m.at(0).at(0) == 1.0);
  CHECK(m.at(0).at(1) == 0.0);
  CHECK(m.at(1).at(0) == 0.0);
  CHECK(m.at(1).at(1) == 1.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::from_2darray()") {
  auto m = lace::square_matrix<double, 2>::from_2darray(
      std::array{std::array{1.0, 2.0}, std::array{3.0, 4.0}});
  CHECK(m.at(0).at(0) == 1.0);
  CHECK(m.at(0).at(1) == 2.0);
  CHECK(m.at(1).at(0) == 3.0);
  CHECK(m.at(1).at(1) == 4.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::to_2darray()") {
  auto m = lace::square_matrix<double, 2>::identity().to_2darray();
  CHECK(m[0][0] == 1.0);
  CHECK(m[0][1] == 0.0);
  CHECK(m[1][0] == 0.0);
  CHECK(m[1][1] == 1.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::to_1darray()") {
  auto m = lace::square_matrix<double, 2>::identity().to_1darray();
  CHECK(m[0] == 1.0);
  CHECK(m[1] == 0.0);
  CHECK(m[2] == 0.0);
  CHECK(m[3] == 1.0);
}

consteval auto matrix_row_nonconst() {
  lace::square_matrix<double, 1> m{1.0};
  m.row(0) = lace::vector<double, 1>{2.0};
  return m[0][0];
}

consteval auto matrix_row_const() {
  const lace::square_matrix<double, 1> m{1.0};
  return m.row(0)[0];
}

TEST_CASE("lace/matrix.hpp: square_matrix::row()") {
  SUBCASE("nonconst") { CHECK(matrix_op_nonconst() == 2.0); }
  SUBCASE("const") { CHECK(matrix_op_const() == 1.0); }
}

consteval auto matrix_col_const() {
  const lace::square_matrix<double, 1> m{1.0};
  return m.col(0)[0];
}

TEST_CASE("lace/matrix.hpp: square_matrix::col()") {
  CHECK(matrix_col_const() == 1.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::transpose()") {
  auto m = lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0}.transpose();
  CHECK(m.at(0).at(0) == 1.0);
  CHECK(m.at(0).at(1) == 3.0);
  CHECK(m.at(1).at(0) == 2.0);
  CHECK(m.at(1).at(1) == 4.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::norm_rows()") {
  auto mn = lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0}.norm_rows();
  CHECK(mn == 7.0);
  mn = lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0}
           .transpose()
           .norm_rows();
  CHECK(mn == 6.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::norm_cols()") {
  auto mn = lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0}.norm_cols();
  CHECK(mn == 6.0);
  mn = lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0}
           .transpose()
           .norm_cols();
  CHECK(mn == 7.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::operator+(vector)") {
  auto m = lace::square_matrix<double, 2>::identity() +
           lace::vector<double, 2>::basis(0);
  CHECK(m.at(0).at(0) == 2.0);
  CHECK(m.at(0).at(1) == 0.0);
  CHECK(m.at(1).at(0) == 1.0);
  CHECK(m.at(1).at(1) == 1.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::operator+(square_matrix)") {
  auto m = lace::square_matrix<double, 2>::identity() +
           lace::square_matrix<double, 2>::identity();
  CHECK(m.at(0).at(0) == 2.0);
  CHECK(m.at(0).at(1) == 0.0);
  CHECK(m.at(1).at(0) == 0.0);
  CHECK(m.at(1).at(1) == 2.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::operator-(vector)") {
  auto m = lace::square_matrix<double, 2>::identity() -
           lace::vector<double, 2>::basis(0);
  CHECK(m.at(0).at(0) == 0.0);
  CHECK(m.at(0).at(1) == 0.0);
  CHECK(m.at(1).at(0) == -1.0);
  CHECK(m.at(1).at(1) == 1.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::operator-(square_matrix)") {
  auto m = lace::square_matrix<double, 2>::identity() -
           lace::square_matrix<double, 2>::identity();
  CHECK(m.at(0).at(0) == 0.0);
  CHECK(m.at(0).at(1) == 0.0);
  CHECK(m.at(1).at(0) == 0.0);
  CHECK(m.at(1).at(1) == 0.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::operator*(double)") {
  auto m = lace::square_matrix<double, 2>::identity() * 2.0;
  CHECK(m.at(0).at(0) == 2.0);
  CHECK(m.at(0).at(1) == 0.0);
  CHECK(m.at(1).at(0) == 0.0);
  CHECK(m.at(1).at(1) == 2.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::operator*(vector)") {
  auto m = lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0} *
           lace::vector<double, 2>::basis(0);
  CHECK(m.at(0) == 1.0);
  CHECK(m.at(1) == 3.0);
  m = lace::vector<double, 2>::basis(0) *
      lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0};
  CHECK(m.at(0) == 1.0);
  CHECK(m.at(1) == 2.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::operator*(square_matrix)") {
  auto m = lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0} *
           lace::square_matrix<double, 2>{4.0, 3.0, 2.0, 1.0};
  CHECK(m.at(0).at(0) == 8.0);
  CHECK(m.at(0).at(1) == 5.0);
  CHECK(m.at(1).at(0) == 20.0);
  CHECK(m.at(1).at(1) == 13.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::swap()") {
  auto m = lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0}.swap(0, 1);
  CHECK(m.at(0).at(0) == 3.0);
  CHECK(m.at(0).at(1) == 4.0);
  CHECK(m.at(1).at(0) == 1.0);
  CHECK(m.at(1).at(1) == 2.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::find_pivot()") {
  auto p = lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0}.find_pivot(0);
  CHECK(p == 1.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::PLU()") {
  auto t = lace::square_matrix<double, 2>{1.0, 2.0, 3.0, 4.0}.PLU();
  auto m = std::get<0>(t);
  CHECK(m.at(0).at(0) == 0.0);
  CHECK(m.at(0).at(1) == 1.0);
  CHECK(m.at(1).at(0) == 1.0);
  CHECK(m.at(1).at(1) == 0.0);
  m = std::get<1>(t);
  CHECK(m.at(0).at(0) - 3.0 == 0.0);
  CHECK(m.at(0).at(1) == 0.0);
  CHECK(m.at(1).at(0) - 1.0 == 0.0);
  auto v = m.at(1).at(1) - 2.0 / 3.0 + 1.0;
  CHECK(v - 1.0 == 0.0);
  m = std::get<2>(t);
  CHECK(m.at(0).at(0) - 1.0 == 0.0);
  CHECK(m.at(0).at(1) - 4.0 / 3.0 == 0.0);
  CHECK(m.at(1).at(0) == 0.0);
  CHECK(m.at(1).at(1) - 1.0 == 0.0);
  m = std::get<3>(t);
  CHECK(m.at(0).at(0) - 1.0 / 3.0 == 0.0);
  CHECK(m.at(0).at(1) == 0.0);
  v = m.at(1).at(0) + 0.5 + 1.0;
  CHECK(v - 1.0 == 0.0);
  v = m.at(1).at(1) - 1.5 + 4.0;
  CHECK(v - 4.0 == 0.0);
  m = std::get<4>(t);
  CHECK(m.at(0).at(0) - 1.0 == 0.0);
  CHECK(m.at(0).at(1) + 4.0 / 3.0 == 0.0);
  CHECK(m.at(1).at(0) == 0.0);
  CHECK(m.at(1).at(1) - 1.0 == 0.0);
}

TEST_CASE("lace/matrix.hpp: square_matrix::invert()") {
  auto m = lace::square_matrix<double, 3>{1.0, 2.0, 3.0, 2.0, 1.0,
                                          2.0, 3.0, 2.0, 1.0}
               .invert();
  auto v = m.at(0).at(0) + 2.0;
  CHECK(v - 2.0 == -0.375);
  v = m.at(0).at(1) + 2.0;
  CHECK(v - 2.0 == 0.5);
  v = m.at(0).at(2) + 2.0;
  CHECK(v - 2.0 == 0.125);
  v = m.at(1).at(0) + 2.0;
  CHECK(v - 2.0 == 0.5);
  v = m.at(1).at(1) + 2.0;
  CHECK(v - 2.0 == -1.0);
  v = m.at(1).at(2) + 2.0;
  CHECK(v - 2.0 == 0.5);
  v = m.at(2).at(0) + 2.0;
  CHECK(v - 2.0 == 0.125);
  v = m.at(2).at(1) + 2.0;
  CHECK(v - 2.0 == 0.5);
  v = m.at(2).at(2) + 2.0;
  CHECK(v - 2.0 == -0.375);
}
