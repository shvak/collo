#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "astro/kepler.hpp"
#include <array>
#include <doctest/doctest.h>
#include <numbers>

using std::numbers::pi;

TEST_CASE("astro/kepler.hpp: anomalies") {
  CHECK(astro::eccentric_anomaly(
            astro::mean_anomaly(astro::eccentric_anomaly_t{1.0}, 0.5), 0.5) +
            2.0 ==
        3.0);
  CHECK(astro::hyperbolic_anomaly(
            astro::mean_anomaly(astro::hyperbolic_anomaly_t{1.0}, 1.5), 1.5) +
            2.0 ==
        3.0);
  CHECK(astro::true_anomaly(
            astro::mean_anomaly(astro::true_anomaly_t{1.0}, 0.5), 0.5) +
            2.0 ==
        3.0);
  CHECK(astro::parabolic_sigma(
            astro::mean_anomaly(astro::parabolic_sigma_t{1.0})) +
            2.0 ==
        3.0);
  CHECK(astro::d_mean_anomaly(astro::eccentric_anomaly_t{pi}, 0.5) + 2.0 ==
        3.5);
  CHECK(astro::d_mean_anomaly(astro::hyperbolic_anomaly_t{0.0}, 1.5) + 2.0 ==
        2.5);
  CHECK(astro::eccentric_anomaly(
            astro::true_anomaly(astro::eccentric_anomaly_t{1.0}, 0.5), 0.5) +
            2.0 ==
        3.0);
  CHECK(astro::hyperbolic_anomaly(
            astro::true_anomaly(astro::hyperbolic_anomaly_t{1.0}, 1.5), 1.5) +
            2.0 ==
        3.0);
}

TEST_CASE("astro/kepler.hpp: true_angle()") {
  CHECK(astro::true_angle(0.0) == 0.0);
  CHECK(astro::true_angle(pi) == pi);
  CHECK(astro::true_angle(2 * pi) == 0.0);
  CHECK(astro::true_angle(-pi) == pi);
  CHECK(astro::true_angle(-2 * pi) == 0.0);
  CHECK(astro::true_angle(pi / 2) == pi / 2);
  CHECK(astro::true_angle(-pi / 2) == 3 * pi / 2);
}

TEST_CASE("astro/kepler.hpp: transformations") {
  double kappa2 = 4 * pi * pi;
  auto ke1 = std::array{1.0 - 0.5 * 0.5, 0.5, 1.0, 1.0, 1.0, 1.0};
  auto ke1_test =
      astro::kepler_elements(astro::cartesian_elements(ke1, kappa2), kappa2);
  CHECK(ke1[0] + 4.0 == ke1_test[0] + 4.0);
  CHECK(ke1[1] + 4.0 == ke1_test[1] + 4.0);
  CHECK(ke1[2] + 4.0 == ke1_test[2] + 4.0);
  CHECK(ke1[3] + 4.0 == ke1_test[3] + 4.0);
  CHECK(ke1[4] + 4.0 == ke1_test[4] + 4.0);
  CHECK(ke1[5] + 4.0 == ke1_test[5] + 4.0);
  auto ke2 = std::array{1.5 * 1.5 - 1.0, 1.5, 1.0, 1.0, 1.0, 1.0};
  auto ke2_test =
      astro::kepler_elements(astro::cartesian_elements(ke2, kappa2), kappa2);
  CHECK(ke2[0] + 4.0 == ke2_test[0] + 4.0);
  CHECK(ke2[1] + 4.0 == ke2_test[1] + 4.0);
  CHECK(ke2[2] + 4.0 == ke2_test[2] + 4.0);
  CHECK(ke2[3] + 4.0 == ke2_test[3] + 4.0);
  CHECK(ke2[4] + 4.0 == ke2_test[4] + 4.0);
  CHECK(ke2[5] + 8.0 == ke2_test[5] + 8.0);
  auto ke3 = std::array{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  auto ke3_test =
      astro::kepler_elements(astro::cartesian_elements(ke3, kappa2), kappa2);
  CHECK(ke3[0] + 4.0 == ke3_test[0] + 4.0);
  CHECK(ke3[1] + 4.0 == ke3_test[1] + 4.0);
  CHECK(ke3[2] + 4.0 == ke3_test[2] + 4.0);
  CHECK(ke3[3] + 4.0 == ke3_test[3] + 4.0);
  CHECK(ke3[4] + 4.0 == ke3_test[4] + 4.0);
  CHECK(ke3[5] + 8.0 == ke3_test[5] + 8.0);
}
