#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "numm/legendre.hpp"
#include <cmath>
#include <doctest/doctest.h>

TEST_CASE("numm/legendre.hpp: legendre<5>(x)") {
  SUBCASE("x=-1") {
    auto t = numm::legendre<5>(-1.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 1.0);
    CHECK(t[1] == -1.0);
    CHECK(t[2] == 1.0);
    CHECK(t[3] == -1.0);
    CHECK(t[4] == 1.0);
    CHECK(t[5] == -1.0);
  }
  SUBCASE("x=0") {
    auto t = numm::legendre<5>(0.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 1.0);
    CHECK(t[1] == 0.0);
    CHECK(t[2] == -0.5);
    CHECK(t[3] == 0.0);
    CHECK(t[4] == 0.375);
    CHECK(t[5] == 0.0);
  }
  SUBCASE("x=1") {
    auto t = numm::legendre<5>(1.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 1.0);
    CHECK(t[1] == 1.0);
    CHECK(t[2] == 1.0);
    CHECK(t[3] == 1.0);
    CHECK(t[4] == 1.0);
    CHECK(t[5] == 1.0);
  }
}

TEST_CASE("numm/legendre.hpp: legendre<5>({-1.0, 0.0, 1.0})") {
  auto t = numm::legendre<5>(std::array<double, 3>{-1.0, 0.0, 1.0});
  CHECK(t.size() == 3);
  CHECK(t[0].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[0][0] == 1.0);
  CHECK(t[0][1] == -1.0);
  CHECK(t[0][2] == 1.0);
  CHECK(t[0][3] == -1.0);
  CHECK(t[0][4] == 1.0);
  CHECK(t[0][5] == -1.0);
  CHECK(t[1][0] == 1.0);
  CHECK(t[1][1] == 0.0);
  CHECK(t[1][2] == -0.5);
  CHECK(t[1][3] == 0.0);
  CHECK(t[1][4] == 0.375);
  CHECK(t[1][5] == 0.0);
  CHECK(t[2][0] == 1.0);
  CHECK(t[2][1] == 1.0);
  CHECK(t[2][2] == 1.0);
  CHECK(t[2][3] == 1.0);
  CHECK(t[2][4] == 1.0);
  CHECK(t[2][5] == 1.0);
}

TEST_CASE("numm/legendre.hpp: dlegendre<5>(x)") {
  SUBCASE("x=-1") {
    auto t = numm::dlegendre<5>(-1.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 1.0);
    CHECK(t[2] == -3.0);
    CHECK(t[3] == 6.0);
    CHECK(t[4] == -10.0);
    CHECK(t[5] == 15.0);
  }
  SUBCASE("x=0") {
    auto t = numm::dlegendre<5>(0.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 1.0);
    CHECK(t[2] == 0.0);
    CHECK(t[3] == -1.5);
    CHECK(t[4] == 0.0);
    CHECK(t[5] == 1.875);
  }
  SUBCASE("x=1") {
    auto t = numm::dlegendre<5>(1.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 1.0);
    CHECK(t[2] == 3.0);
    CHECK(t[3] == 6.0);
    CHECK(t[4] == 10.0);
    CHECK(t[5] == 15.0);
  }
}

TEST_CASE("numm/legendre.hpp: dlegendre<5>({-1.0, 0.0, 1.0})") {
  auto t = numm::dlegendre<5>(std::array<double, 3>{-1.0, 0.0, 1.0});
  CHECK(t.size() == 3);
  CHECK(t[0].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[0][0] == 0.0);
  CHECK(t[0][1] == 1.0);
  CHECK(t[0][2] == -3.0);
  CHECK(t[0][3] == 6.0);
  CHECK(t[0][4] == -10.0);
  CHECK(t[0][5] == 15.0);
  CHECK(t[1][0] == 0.0);
  CHECK(t[1][1] == 1.0);
  CHECK(t[1][2] == 0.0);
  CHECK(t[1][3] == -1.5);
  CHECK(t[1][4] == 0.0);
  CHECK(t[1][5] == 1.875);
  CHECK(t[2][0] == 0.0);
  CHECK(t[2][1] == 1.0);
  CHECK(t[2][2] == 3.0);
  CHECK(t[2][3] == 6.0);
  CHECK(t[2][4] == 10.0);
  CHECK(t[2][5] == 15.0);
}

TEST_CASE("numm/legendre.hpp: d2legendre<5>(x)") {
  SUBCASE("x=-1") {
    auto t = numm::d2legendre<5>(-1.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 0.0);
    CHECK(t[2] == 3.0);
    CHECK(t[3] == -15.0);
    CHECK(t[4] == 45.0);
    CHECK(t[5] == -105.0);
  }
  SUBCASE("x=0") {
    auto t = numm::d2legendre<5>(0.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 0.0);
    CHECK(t[2] == 3.0);
    CHECK(t[3] == 0.0);
    CHECK(t[4] == -7.5);
    CHECK(t[5] == 0.0);
  }
  SUBCASE("x=1") {
    auto t = numm::d2legendre<5>(1.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 0.0);
    CHECK(t[2] == 3.0);
    CHECK(t[3] == 15.0);
    CHECK(t[4] == 45.0);
    CHECK(t[5] == 105.0);
  }
}

TEST_CASE("numm/legendre.hpp: d2legendre<5>({-1.0, 0.0, 1.0})") {
  auto t = numm::d2legendre<5>(std::array<double, 3>{-1.0, 0.0, 1.0});
  CHECK(t.size() == 3);
  CHECK(t[0].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[0][0] == 0.0);
  CHECK(t[0][1] == 0.0);
  CHECK(t[0][2] == 3.0);
  CHECK(t[0][3] == -15.0);
  CHECK(t[0][4] == 45.0);
  CHECK(t[0][5] == -105.0);
  CHECK(t[1][0] == 0.0);
  CHECK(t[1][1] == 0.0);
  CHECK(t[1][2] == 3.0);
  CHECK(t[1][3] == 0.0);
  CHECK(t[1][4] == -7.5);
  CHECK(t[1][5] == 0.0);
  CHECK(t[2][0] == 0.0);
  CHECK(t[2][1] == 0.0);
  CHECK(t[2][2] == 3.0);
  CHECK(t[2][3] == 15.0);
  CHECK(t[2][4] == 45.0);
  CHECK(t[2][5] == 105.0);
}

TEST_CASE("numm/legendre.hpp: legendre_sh<5>(x)") {
  SUBCASE("x=0") {
    auto t = numm::legendre_sh<5>(0.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 1.0);
    CHECK(t[1] == -1.0);
    CHECK(t[2] == 1.0);
    CHECK(t[3] == -1.0);
    CHECK(t[4] == 1.0);
    CHECK(t[5] == -1.0);
  }
  SUBCASE("x=0.5") {
    auto t = numm::legendre_sh<5>(0.5);
    CHECK(t.size() == 6);
    CHECK(t[0] == 1.0);
    CHECK(t[1] == 0.0);
    CHECK(t[2] == -0.5);
    CHECK(t[3] == 0.0);
    CHECK(t[4] == 0.375);
    CHECK(t[5] == 0.0);
  }
  SUBCASE("x=1") {
    auto t = numm::legendre_sh<5>(1.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 1.0);
    CHECK(t[1] == 1.0);
    CHECK(t[2] == 1.0);
    CHECK(t[3] == 1.0);
    CHECK(t[4] == 1.0);
    CHECK(t[5] == 1.0);
  }
}

TEST_CASE("numm/legendre.hpp: legendre_sh<5>({0.0, 0.5, 1.0})") {
  auto t = numm::legendre_sh<5>(std::array<double, 3>{0.0, 0.5, 1.0});
  CHECK(t.size() == 3);
  CHECK(t[0].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[0][0] == 1.0);
  CHECK(t[0][1] == -1.0);
  CHECK(t[0][2] == 1.0);
  CHECK(t[0][3] == -1.0);
  CHECK(t[0][4] == 1.0);
  CHECK(t[0][5] == -1.0);
  CHECK(t[1][0] == 1.0);
  CHECK(t[1][1] == 0.0);
  CHECK(t[1][2] == -0.5);
  CHECK(t[1][3] == 0.0);
  CHECK(t[1][4] == 0.375);
  CHECK(t[1][5] == 0.0);
  CHECK(t[2][0] == 1.0);
  CHECK(t[2][1] == 1.0);
  CHECK(t[2][2] == 1.0);
  CHECK(t[2][3] == 1.0);
  CHECK(t[2][4] == 1.0);
  CHECK(t[2][5] == 1.0);
}

TEST_CASE("numm/legendre.hpp: dlegendre_sh<5>(x)") {
  SUBCASE("x=-1") {
    auto t = numm::dlegendre_sh<5>(0.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 2.0);
    CHECK(t[2] == -6.0);
    CHECK(t[3] == 12.0);
    CHECK(t[4] == -20.0);
    CHECK(t[5] == 30.0);
  }
  SUBCASE("x=0") {
    auto t = numm::dlegendre_sh<5>(0.5);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 2.0);
    CHECK(t[2] == 0.0);
    CHECK(t[3] == -3.0);
    CHECK(t[4] == 0.0);
    CHECK(t[5] == 3.75);
  }
  SUBCASE("x=1") {
    auto t = numm::dlegendre_sh<5>(1.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 2.0);
    CHECK(t[2] == 6.0);
    CHECK(t[3] == 12.0);
    CHECK(t[4] == 20.0);
    CHECK(t[5] == 30.0);
  }
}

TEST_CASE("numm/legendre.hpp: dlegendre_sh<5>({0.0, 0.5, 1.0})") {
  auto t = numm::dlegendre_sh<5>(std::array<double, 3>{0.0, 0.5, 1.0});
  CHECK(t.size() == 3);
  CHECK(t[0].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[0][0] == 0.0);
  CHECK(t[0][1] == 2.0);
  CHECK(t[0][2] == -6.0);
  CHECK(t[0][3] == 12.0);
  CHECK(t[0][4] == -20.0);
  CHECK(t[0][5] == 30.0);
  CHECK(t[1][0] == 0.0);
  CHECK(t[1][1] == 2.0);
  CHECK(t[1][2] == 0.0);
  CHECK(t[1][3] == -3.0);
  CHECK(t[1][4] == 0.0);
  CHECK(t[1][5] == 3.75);
  CHECK(t[2][0] == 0.0);
  CHECK(t[2][1] == 2.0);
  CHECK(t[2][2] == 6.0);
  CHECK(t[2][3] == 12.0);
  CHECK(t[2][4] == 20.0);
  CHECK(t[2][5] == 30.0);
}

TEST_CASE("numm/legendre.hpp: d2legendre_sh<5>(x)") {
  SUBCASE("x=-1") {
    auto t = numm::d2legendre_sh<5>(0.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 0.0);
    CHECK(t[2] == 12.0);
    CHECK(t[3] == -60.0);
    CHECK(t[4] == 180.0);
    CHECK(t[5] == -420.0);
  }
  SUBCASE("x=0") {
    auto t = numm::d2legendre_sh<5>(0.5);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 0.0);
    CHECK(t[2] == 12.0);
    CHECK(t[3] == 0.0);
    CHECK(t[4] == -30.0);
    CHECK(t[5] == 0.0);
  }
  SUBCASE("x=1") {
    auto t = numm::d2legendre_sh<5>(1.0);
    CHECK(t.size() == 6);
    CHECK(t[0] == 0.0);
    CHECK(t[1] == 0.0);
    CHECK(t[2] == 12.0);
    CHECK(t[3] == 60.0);
    CHECK(t[4] == 180.0);
    CHECK(t[5] == 420.0);
  }
}

TEST_CASE("numm/legendre.hpp: d2legendre_sh<5>({0.0, 0.5, 1.0})") {
  auto t = numm::d2legendre_sh<5>(std::array<double, 3>{0.0, 0.5, 1.0});
  CHECK(t.size() == 3);
  CHECK(t[0].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[1].size() == 6);
  CHECK(t[0][0] == 0.0);
  CHECK(t[0][1] == 0.0);
  CHECK(t[0][2] == 12.0);
  CHECK(t[0][3] == -60.0);
  CHECK(t[0][4] == 180.0);
  CHECK(t[0][5] == -420.0);
  CHECK(t[1][0] == 0.0);
  CHECK(t[1][1] == 0.0);
  CHECK(t[1][2] == 12.0);
  CHECK(t[1][3] == 0.0);
  CHECK(t[1][4] == -30.0);
  CHECK(t[1][5] == 0.0);
  CHECK(t[2][0] == 0.0);
  CHECK(t[2][1] == 0.0);
  CHECK(t[2][2] == 12.0);
  CHECK(t[2][3] == 60.0);
  CHECK(t[2][4] == 180.0);
  CHECK(t[2][5] == 420.0);
}

TEST_CASE("numm/legendre.hpp: roots_legendre<3, double>()") {
  auto t = numm::roots_legendre<3, double>();
  CHECK(t.size() == 3);
  auto r = std::sqrt(0.6);
  t[0] += r + 1.0;
  t[1] += 1.0;
  t[2] -= r - 1.0;
  CHECK(t[0] - 1.0 == 0.0);
  CHECK(t[1] - 1.0 == 0.0);
  CHECK(t[2] - 1.0 == 0.0);
}

TEST_CASE("numm/legendre.hpp: roots_dlegendre<3, double>()") {
  auto t = numm::roots_dlegendre<3, double>();
  CHECK(t.size() == 2);
  auto r = std::sqrt(0.2);
  t[0] += r + 1.0;
  t[1] -= r - 1.0;
  CHECK(t[0] - 1.0 == 0.0);
  CHECK(t[1] - 1.0 == 0.0);
}

TEST_CASE("numm/legendre.hpp: roots_ilegendre<3, double>()") {
  auto t = numm::roots_ilegendre<3, double>();
  CHECK(t.size() == 4);
  auto r = std::sqrt(0.2);
  t[1] += r + 1.0;
  t[2] -= r - 1.0;
  CHECK(t[0] == -1.0);
  CHECK(t[1] - 1.0 == 0.0);
  CHECK(t[2] - 1.0 == 0.0);
  CHECK(t[3] == 1.0);
}

TEST_CASE("numm/legendre.hpp: roots_legendre_sh<3, double>()") {
  auto t = numm::roots_legendre_sh<3, double>();
  CHECK(t.size() == 3);
  auto r = std::sqrt(0.6);
  t[0] += r / 2 + 1.0;
  t[1] += 1.0;
  t[2] -= r / 2 - 1.0;
  CHECK(t[0] - 1.0 == 0.5);
  CHECK(t[1] - 1.0 == 0.5);
  CHECK(t[2] - 1.0 == 0.5);
}

TEST_CASE("numm/legendre.hpp: roots_dlegendre_sh<3, double>()") {
  auto t = numm::roots_dlegendre_sh<3, double>();
  CHECK(t.size() == 2);
  auto r = std::sqrt(0.2);
  t[0] += r / 2 + 1.0;
  t[1] -= r / 2 - 1.0;
  CHECK(t[0] - 1.0 == 0.5);
  CHECK(t[1] - 1.0 == 0.5);
}

TEST_CASE("numm/legendre.hpp: roots_ilegendre_sh<3, double>()") {
  auto t = numm::roots_ilegendre_sh<3, double>();
  CHECK(t.size() == 4);
  auto r = std::sqrt(0.2);
  t[1] += r / 2 + 1.0;
  t[2] -= r / 2 - 1.0;
  CHECK(t[0] == 0.0);
  CHECK(t[1] - 1.0 == 0.5);
  CHECK(t[2] - 1.0 == 0.5);
  CHECK(t[3] == 1.0);
}
