#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "numm/base.hpp"
#include <array>
#include <cmath>
#include <doctest/doctest.h>
#include <numbers>

TEST_CASE("numm/base.hpp: cos(x)") {
  auto pi = std::numbers::pi_v<double>;

  SUBCASE("x=0") {
    CHECK(numm::cos(0.0) == doctest::Approx(1.0).epsilon(1e-15));
  }
  SUBCASE("x=pi") {
    CHECK(numm::cos(pi) == doctest::Approx(-1.0).epsilon(1e-15));
  }
  SUBCASE("x=-pi") {
    CHECK(numm::cos(-pi) == doctest::Approx(-1.0).epsilon(1e-15));
  }
  SUBCASE("x=pi/2") {
    CHECK(numm::cos(pi / 2) + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  }
  SUBCASE("x=-pi/2") {
    CHECK(numm::cos(-pi / 2) + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  }
}

// accuracy of numm::sqrt() isn't the best
TEST_CASE("numm/base.hpp: sqrt(x)") {
  SUBCASE("x=1") {
    CHECK(numm::sqrt(1.0) == doctest::Approx(1.0).epsilon(1e-15));
  }
  SUBCASE("x=2") {
    CHECK(numm::sqrt(2.0) == doctest::Approx(std::sqrt(2.0)).epsilon(1e-15));
  }
  SUBCASE("x=4") {
    CHECK(numm::sqrt(4.0) == doctest::Approx(2.0).epsilon(1e-15));
  }
}

TEST_CASE("numm/base.hpp: cutoff_first<n>({1, 2, 3})") {
  std::array<double, 3> arr = {1.0, 2.0, 3.0};
  auto res0 = numm::cutoff_first<0>(arr);
  auto res1 = numm::cutoff_first<1>(arr);
  auto res2 = numm::cutoff_first<2>(arr);
  SUBCASE("n=0") {
    CHECK(res0.size() == 3);
    CHECK(res0[0] == 1.0);
    CHECK(res0[1] == 2.0);
    CHECK(res0[2] == 3.0);
  }
  SUBCASE("n=1") {
    CHECK(res1.size() == 2);
    CHECK(res1[0] == 2.0);
    CHECK(res1[1] == 3.0);
  }
  SUBCASE("n=2") {
    CHECK(res2.size() == 1);
    CHECK(res2[0] == 3.0);
  }
}

TEST_CASE("numm/base.hpp: ortho_poly<5, 0, double>(...)") {
  auto res =
      numm::ortho_poly<5, 0, double>(1.0, {0.0, 1.0}, [](double, auto &arr) {
        for (std::size_t i = 1; i < 5; ++i)
          arr[i + 1] = arr[i] + arr[i - 1];
      });
  CHECK(res.size() == 6);
  CHECK(res[0] == 0.0);
  CHECK(res[1] == 1.0);
  CHECK(res[2] == 1.0);
  CHECK(res[3] == 2.0);
  CHECK(res[4] == 3.0);
  CHECK(res[5] == 5.0);
}

TEST_CASE("numm/base.hpp: ortho_poly<3, 1, 3, double>(...)") {
  auto res = numm::ortho_poly<3, 1, 3, double>({0.0, 1.0, 2.0}, [](double) {
    return std::array<double, 3>{1.0, 2.0, 3.0};
  });
  CHECK(res.size() == 3);
  CHECK(res[0].size() == 3);
  CHECK(res[0][0] == 1.0);
  CHECK(res[0][1] == 2.0);
  CHECK(res[0][2] == 3.0);
  CHECK(res[1].size() == 3);
  CHECK(res[1][0] == 1.0);
  CHECK(res[1][1] == 2.0);
  CHECK(res[1][2] == 3.0);
  CHECK(res[2].size() == 3);
  CHECK(res[2][0] == 1.0);
  CHECK(res[2][1] == 2.0);
  CHECK(res[2][2] == 3.0);
}
