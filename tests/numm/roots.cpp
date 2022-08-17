#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "numm/roots.hpp"
#include <doctest/doctest.h>

double f1(double x) { return 0.25 - x * x; }

double f2(double x) { return 1.00 - x * x; }

double f3(double x) { return 4.00 - x * x; }

double df(double x) { return -2 * x; }

TEST_CASE("numm/roots.hpp: newton(...)") {
  SUBCASE("f=a^2-x^2") {
    CHECK(numm::newton(1.0, f1, df) == 0.5);
    CHECK(numm::newton(2.0, f2, df) == 1.0);
    CHECK(numm::newton(3.0, f3, df) == 2.0);
  }
  SUBCASE("f=x^2-a^2") {
    auto g1 = [](double x) { return x * x - 0.25; };
    auto g2 = [](double x) { return x * x - 1.00; };
    auto g3 = [](double x) { return x * x - 4.00; };
    auto dg = [](double x) { return 2.0 * x; };
    CHECK(numm::newton(1.0, g1, dg) == 0.5);
    CHECK(numm::newton(2.0, g2, dg) == 1.0);
    CHECK(numm::newton(3.0, g3, dg) == 2.0);
  }
}

TEST_CASE("numm/roots.hpp: fixed_point(...)") {
  SUBCASE("f=a^2-x^2") {
    auto t = numm::fixed_point(1.0, f1) + 1.0;
    CHECK(t - 1.0 == 0.5);
    t = numm::fixed_point(2.0, f2) + 1.0;
    CHECK(t - 1.0 == 1.0);
    t = numm::fixed_point(3.0, f3) + 1.0;
    CHECK(t - 1.0 == 2.0);
  }
  SUBCASE("f=x^2-a^2") {
    auto g1 = [](double x) { return x * x - 0.25; };
    auto g2 = [](double x) { return x * x - 1.00; };
    auto g3 = [](double x) { return x * x - 4.00; };
    auto t = numm::fixed_point(1.0, g1) + 1.0;
    CHECK(t - 1.0 == 0.5);
    t = numm::fixed_point(2.0, g2) + 1.0;
    CHECK(t - 1.0 == 1.0);
    t = numm::fixed_point(3.0, g3) + 1.0;
    CHECK(t - 1.0 == 2.0);
  }
}
