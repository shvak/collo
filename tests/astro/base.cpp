#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "astro/base.hpp"
#include <doctest/doctest.h>

using namespace std::chrono;

TEST_CASE("astro/base.hpp: julian_date()") {
  CHECK(astro::julian_date<double>(29d / 01 / 1983, 0h + 0min + 0s) ==
        2445363.5);
}

TEST_CASE("astro/base.hpp: gregorian_from_jd()") {
  for (int t = 0; t < 10; ++t) {
    CAPTURE(t);
    auto ymd = astro::gregorian_from_jd(2445363.5 + t * .1);
    CHECK(static_cast<int>(ymd.year()) == 1983);
    CHECK(static_cast<unsigned>(ymd.month()) == 1);
    CHECK(static_cast<unsigned>(ymd.day()) == 29);
  }
}

TEST_CASE("astro/base.hpp: time_from_jd()") {
  for (int t = 0; t < 10; ++t) {
    CAPTURE(t);
    auto time = astro::time_from_jd(2445363.5 + t * .1);
    CHECK(static_cast<int>(time.count() * 10 + .1) == t);
  }
}
