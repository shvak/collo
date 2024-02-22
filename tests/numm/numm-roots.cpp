#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "numm/roots.hpp"
#include <doctest/doctest.h>

double f1(double x) { return 0.25 - x * x; }

double f2(double x) { return 1.00 - x * x; }

double f3(double x) { return 4.00 - x * x; }

double df(double x) { return -2 * x; }

double g1(double x) { return x * x - 0.25; }

double g2(double x) { return x * x - 1.0; }

double g3(double x) { return x * x - 4.0; }

double dg(double x) { return 2 * x; }

TEST_CASE("numm/roots.hpp: newton(...)") {
  SUBCASE("f=a^2-x^2") {
    CHECK(numm::newton(1.0, f1, df) == doctest::Approx(0.5).epsilon(1e-15));
    CHECK(numm::newton(2.0, f2, df) == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(numm::newton(3.0, f3, df) == doctest::Approx(2.0).epsilon(1e-15));
  }
  SUBCASE("f=x^2-a^2") {
    CHECK(numm::newton(1.0, g1, dg) == doctest::Approx(0.5).epsilon(1e-15));
    CHECK(numm::newton(2.0, g2, dg) == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(numm::newton(3.0, g3, dg) == doctest::Approx(2.0).epsilon(1e-15));
  }
}

TEST_CASE("numm/roots.hpp: fixed_point(...)") {
  SUBCASE("f=a^2-x^2") {
    CHECK(numm::fixed_point(1.0, f1) == doctest::Approx(0.5).epsilon(1e-15));
    CHECK(numm::fixed_point(2.0, f2) == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(numm::fixed_point(3.0, f3) == doctest::Approx(2.0).epsilon(1e-15));
  }
  SUBCASE("f=x^2-a^2") {
    CHECK(numm::fixed_point(1.0, g1) == doctest::Approx(0.5).epsilon(1e-15));
    CHECK(numm::fixed_point(2.0, g2) == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(numm::fixed_point(3.0, g3) == doctest::Approx(2.0).epsilon(1e-15));
  }
}
