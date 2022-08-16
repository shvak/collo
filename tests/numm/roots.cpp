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
    auto f1 = [](double x) { return x * x - 0.25; };
    auto f2 = [](double x) { return x * x - 1.00; };
    auto f3 = [](double x) { return x * x - 4.00; };
    auto df = [](double x) { return 2.0 * x; };
    CHECK(numm::newton(1.0, f1, df) == 0.5);
    CHECK(numm::newton(2.0, f2, df) == 1.0);
    CHECK(numm::newton(3.0, f3, df) == 2.0);
  }
}
