#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "astro/kepler.hpp"
#include <array>
#include <doctest/doctest.h>
#include <numbers>

using std::numbers::pi;

TEST_CASE("astro/kepler.hpp: anomalies") {
  CHECK(astro::eccentric_anomaly(
            astro::mean_anomaly(astro::eccentric_anomaly_t<double>{1.0}, 0.5),
            0.5) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::hyperbolic_anomaly(
            astro::mean_anomaly(astro::hyperbolic_anomaly_t<double>{1.0}, 1.5),
            1.5) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::true_anomaly(
            astro::mean_anomaly(astro::true_anomaly_t<double>{1.0}, 0.5),
            0.5) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::parabolic_sigma(
            astro::mean_anomaly(astro::parabolic_sigma_t<double>{1.0})) ==
        doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::d_mean_anomaly(astro::eccentric_anomaly_t<double>{pi}, 0.5) ==
        doctest::Approx(1.5).epsilon(1e-15));
  CHECK(astro::d_mean_anomaly(astro::hyperbolic_anomaly_t<double>{0.0}, 1.5) ==
        doctest::Approx(0.5).epsilon(1e-15));
  CHECK(astro::eccentric_anomaly(
            astro::true_anomaly(astro::eccentric_anomaly_t<double>{1.0}, 0.5),
            0.5) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::hyperbolic_anomaly(
            astro::true_anomaly(astro::hyperbolic_anomaly_t<double>{1.0}, 1.5),
            1.5) == doctest::Approx(1.0).epsilon(1e-15));
}

TEST_CASE("astro/kepler.hpp: true_angle()") {
  CHECK(astro::true_angle(0.0) + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::true_angle(pi) - pi + 1.0 ==
        doctest::Approx(1.0).epsilon(1e-15));
  // CHECK(astro::true_angle(2 * pi) + 1.0 ==
  // doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::true_angle(-pi) - pi + 1.0 ==
        doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::true_angle(-2 * pi) + 1.0 ==
        doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::true_angle(pi / 2) - pi / 2 + 1.0 ==
        doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::true_angle(-pi / 2) - 3 * pi / 2 + 1.0 ==
        doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::true_angle(3 * pi / 2) - 3 * pi / 2 + 1.0 ==
        doctest::Approx(1.0).epsilon(1e-15));
  CHECK(astro::true_angle(-3 * pi / 2) - pi / 2 + 1.0 ==
        doctest::Approx(1.0).epsilon(1e-15));
}

TEST_CASE("astro/kepler.hpp: transformations") {
  double kappa2 = 4 * pi * pi;
  auto ke1 = std::array{1.0 - 0.5 * 0.5, 0.5, 1.0, 1.0, 1.0, 1.0};
  auto ke1_test =
      astro::kepler_elements(astro::cartesian_elements(ke1, kappa2), kappa2);
  CHECK(ke1[0] - ke1_test[0] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke1[1] - ke1_test[1] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke1[2] - ke1_test[2] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke1[3] - ke1_test[3] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke1[4] - ke1_test[4] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke1[5] - ke1_test[5] + 1.0 == doctest::Approx(1.0).epsilon(1e-14));
  auto ke2 = std::array{1.5 * 1.5 - 1.0, 1.5, 1.0, 1.0, 1.0, 1.0};
  auto ke2_test =
      astro::kepler_elements(astro::cartesian_elements(ke2, kappa2), kappa2);
  CHECK(ke2[0] - ke2_test[0] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke2[1] - ke2_test[1] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke2[2] - ke2_test[2] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke2[3] - ke2_test[3] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke2[4] - ke2_test[4] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke2[5] - ke2_test[5] + 1.0 == doctest::Approx(1.0).epsilon(1e-14));
  auto ke3 = std::array{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  auto ke3_test =
      astro::kepler_elements(astro::cartesian_elements(ke3, kappa2), kappa2);
  CHECK(ke3[0] - ke3_test[0] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke3[1] - ke3_test[1] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke3[2] - ke3_test[2] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke3[3] - ke3_test[3] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke3[4] - ke3_test[4] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(ke3[5] - ke3_test[5] + 1.0 == doctest::Approx(1.0).epsilon(1e-14));
  for (double M = 0.; M < 2 * pi; M += .1) {
    auto ke4 = std::array{0.2, 0.9, 1.0, 1.0, 1.0, M};
    auto ke4_test =
        astro::kepler_elements(astro::cartesian_elements(ke4, kappa2), kappa2);
    CAPTURE(M);
    CHECK(ke4[0] - ke4_test[0] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(ke4[1] - ke4_test[1] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(ke4[2] - ke4_test[2] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(ke4[3] - ke4_test[3] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(ke4[4] - ke4_test[4] + 1.0 == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(ke4[5] - ke4_test[5] + 1.0 == doctest::Approx(1.0).epsilon(1e-14));
  }
}
