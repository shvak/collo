#ifndef ASTRO_KEPLER_HPP
#define ASTRO_KEPLER_HPP

#include "base.hpp"
#include <cmath>
#include <numbers>
#include <numm/roots.hpp>

namespace astro {

template <numm::number_type num_t>
constexpr mean_anomaly_t<num_t> mean_anomaly(eccentric_anomaly_t<num_t> E,
                                             num_t e) {
  return E - e * std::sin(E);
}

template <numm::number_type num_t>
constexpr mean_anomaly_t<num_t> mean_anomaly(hyperbolic_anomaly_t<num_t> H,
                                             num_t e) {
  return e * std::sinh(H) - H;
}

template <numm::number_type num_t>
constexpr mean_anomaly_t<num_t> mean_anomaly(true_anomaly_t<num_t> theta,
                                             num_t e) {
  return mean_anomaly(eccentric_anomaly(theta, e), e);
}

template <numm::number_type num_t>
constexpr mean_anomaly_t<num_t> mean_anomaly(parabolic_sigma_t<num_t> sigma) {
  return (sigma * sigma + num_t{3.0}) * sigma / num_t{2.0};
}

template <numm::number_type num_t>
constexpr d_mean_anomaly_t<num_t> d_mean_anomaly(eccentric_anomaly_t<num_t> E,
                                                 num_t e) {
  return num_t{1.0} - e * std::cos(E);
}

template <numm::number_type num_t>
constexpr d_mean_anomaly_t<num_t> d_mean_anomaly(hyperbolic_anomaly_t<num_t> H,
                                                 num_t e) {
  return e * std::cosh(H) - num_t{1.0};
}

template <numm::number_type num_t>
constexpr eccentric_anomaly_t<num_t> eccentric_anomaly(mean_anomaly_t<num_t> M,
                                                       num_t e) {
  auto f = [M, e](num_t x) {
    return mean_anomaly(eccentric_anomaly_t{x}, e) - M;
  };
  auto df = [e](num_t x) { return d_mean_anomaly(eccentric_anomaly_t{x}, e); };
  return numm::newton(static_cast<num_t>(M), f, df);
}

template <numm::number_type num_t>
constexpr eccentric_anomaly_t<num_t>
eccentric_anomaly(true_anomaly_t<num_t> theta, num_t e) {
  num_t beta = e / (1 + std::sqrt(1 - e * e));
  return theta -
         2 * std::atan2(beta * std::sin(theta), 1 + beta * std::cos(theta));
}

template <numm::number_type num_t>
constexpr hyperbolic_anomaly_t<num_t>
hyperbolic_anomaly(mean_anomaly_t<num_t> M, num_t e) {
  auto f = [M, e](num_t x) {
    return mean_anomaly(hyperbolic_anomaly_t{x}, e) - M;
  };
  auto df = [e](num_t x) { return d_mean_anomaly(hyperbolic_anomaly_t{x}, e); };
  return numm::newton(static_cast<num_t>(M), f, df);
}

template <numm::number_type num_t>
constexpr hyperbolic_anomaly_t<num_t>
hyperbolic_anomaly(true_anomaly_t<num_t> theta, num_t e) {
  num_t c_ = std::sqrt((e - 1) / (e + 1));
  return num_t{2.0} * std::atanh(c_ * std::tan(theta / num_t{2.0}));
}

template <numm::number_type num_t>
constexpr true_anomaly_t<num_t> true_anomaly(eccentric_anomaly_t<num_t> E,
                                             num_t e) {
  num_t beta = e / (1 + std::sqrt(1 - e * e));
  return E + 2 * std::atan2(beta * std::sin(E), 1 - beta * std::cos(E));
}

template <numm::number_type num_t>
constexpr true_anomaly_t<num_t> true_anomaly(hyperbolic_anomaly_t<num_t> H,
                                             num_t e) {
  num_t eta = std::sqrt(e * e - num_t{1.0});
  return std::atan2(eta * std::sinh(H), e - std::cosh(H));
}

template <numm::number_type num_t>
constexpr true_anomaly_t<num_t> true_anomaly(mean_anomaly_t<num_t> M, num_t e) {
  return true_anomaly(eccentric_anomaly(M, e), e);
}

template <numm::number_type num_t>
constexpr parabolic_sigma_t<num_t> parabolic_sigma(mean_anomaly_t<num_t> M) {
  num_t x = std::sqrt(num_t{1.0} + M * M);
  return std::cbrt(x + M) - std::cbrt(x - M);
}

template <numm::number_type num_t, typename arr_t>
constexpr auto unknown_res() {
  auto x = num_t{std::numeric_limits<double>::quiet_NaN()};
  return arr_t{x, x, x, x, x, x};
}

template <numm::number_type num_t> constexpr auto true_angle(num_t angle) {
  constexpr num_t pi = std::numbers::pi_v<num_t>;
  return std::fmod(angle + 2 * pi, 2 * pi);
}

template <numm::number_type num_t, typename arr_t>
constexpr auto
kepler_elements(const arr_t &xyz,
                num_t kappa2) // xyz: x, y, z, dot_x, dot_y, dot_z
{
  num_t x = xyz[0], y = xyz[1], z = xyz[2], dot_x = xyz[3], dot_y = xyz[4],
        dot_z = xyz[5];
  num_t c_1 = y * dot_z - z * dot_y, c_2 = z * dot_x - x * dot_z,
        c_3 = x * dot_y - y * dot_x, c2 = c_1 * c_1 + c_2 * c_2 + c_3 * c_3,
        c = std::sqrt(c2), //v2 = dot_x * dot_x + dot_y * dot_y + dot_z * dot_z,
        r = std::sqrt(x * x + y * y + z * z),
        dot_r = (x * dot_x + y * dot_y + z * dot_z) / r, p = c2 / kappa2,
        est = dot_r * std::sqrt(p / kappa2), est2 = est * est,
        ect = p / r - num_t{1.0}, ect2 = ect * ect;

  num_t //h = v2 / num_t{2.0} - kappa2 / r,
        // a = - kappa2 / h / num_t{2.0}, //num_t{1.0} / (num_t{2.0} / r - v2 /
        // kappa2),
      e = std::sqrt(est2 + ect2), theta = true_angle(std::atan2(est, ect)),
        Omega = true_angle(std::atan2(c_1, -c_2)), i = std::acos(c_3 / c),
        g = true_angle(std::atan2(z, std::sin(i) * (x * std::cos(Omega) +
                                                    y * std::sin(Omega))) -
                       theta);

  num_t M;
  if (e + num_t{8.0} > num_t{9.0})
    M = mean_anomaly(hyperbolic_anomaly(true_anomaly_t{theta}, e), e);
  else if (e + num_t{8.0} < num_t{9.0})
    M = mean_anomaly(eccentric_anomaly(true_anomaly_t{theta}, e), e);
  else
    M = mean_anomaly(parabolic_sigma_t<num_t>{std::tan(theta / num_t{2.0})});

  return arr_t{p, e, i, Omega, g, M};
}

template <numm::number_type num_t, typename arr_t>
constexpr auto rotate_from_orbit_plane(const arr_t &ceop, num_t i, num_t Omega,
                                       num_t g) {
  num_t cosi = std::cos(i), sini = std::sin(i), cosO = std::cos(Omega),
        sinO = std::sin(Omega), cosg = std::cos(g), sing = std::sin(g),
        qxx = cosO * cosg - cosi * sinO * sing,
        qxy = -cosO * sing - cosi * sinO * cosg,
        qyx = sinO * cosg + cosi * cosO * sing,
        qyy = -sinO * sing + cosi * cosO * cosg, qzx = sini * sing,
        qzy = sini * cosg;

  num_t X = ceop[0], Y = ceop[1], dot_X = ceop[3], dot_Y = ceop[4];

  num_t x = qxx * X + qxy * Y, y = qyx * X + qyy * Y, z = qzx * X + qzy * Y,
        dot_x = qxx * dot_X + qxy * dot_Y, dot_y = qyx * dot_X + qyy * dot_Y,
        dot_z = qzx * dot_X + qzy * dot_Y;

  return arr_t{x, y, z, dot_x, dot_y, dot_z};
}

template <numm::number_type num_t, typename arr_t>
constexpr auto cartesian_elements_ellipse_orbit_plane(
    num_t p, num_t e, eccentric_anomaly_t<num_t> E, num_t kappa2) {
  num_t cosE = std::cos(E), sinE = std::sin(E), eta2 = num_t{1.0} - e * e,
        eta = std::sqrt(eta2), a = p / eta2, n = std::sqrt(kappa2 / a / a / a),
        dot_E = n / (num_t{1.0} - e * cosE);

  num_t X = a * (cosE - e), Y = p / eta * sinE, dot_X = -a * sinE * dot_E,
        dot_Y = p / eta * cosE * dot_E;

  return arr_t{X, Y, num_t{0.0}, dot_X, dot_Y, num_t{0.0}};
}

template <numm::number_type num_t, typename arr_t>
constexpr auto
cartesian_elements_ellipse(const arr_t &ke,
                           num_t kappa2) // ke: p, e, i, Omega, g, M
{
  num_t p = ke[0], e = ke[1], i = ke[2], Omega = ke[3], g = ke[4], M = ke[5];
  auto E = eccentric_anomaly(mean_anomaly_t{M}, e);

  return rotate_from_orbit_plane(
      cartesian_elements_ellipse_orbit_plane<num_t, arr_t>(p, e, E, kappa2), i,
      Omega, g);
}

template <numm::number_type num_t, typename arr_t>
constexpr auto cartesian_elements_hyperbola_orbit_plane(
    num_t p, num_t e, hyperbolic_anomaly_t<num_t> H, num_t kappa2) {
  num_t chH = std::cosh(H), shH = std::sinh(H), eta2 = e * e - num_t{1.0},
        eta = std::sqrt(eta2), a = -p / eta2,
        n = std::sqrt(-kappa2 / a / a / a), dot_H = n / (e * chH - num_t{1.0});

  num_t X = a * (chH - e), Y = p / eta * shH, dot_X = a * shH * dot_H,
        dot_Y = p / eta * chH * dot_H;

  return arr_t{X, Y, num_t{0.0}, dot_X, dot_Y, num_t{0.0}};
}

template <numm::number_type num_t, typename arr_t>
constexpr auto
cartesian_elements_hyperbola(const arr_t &ke,
                             num_t kappa2) // ke: p, e, i, Omega, g, M
{
  num_t p = ke[0], e = ke[1], i = ke[2], Omega = ke[3], g = ke[4], M = ke[5];
  auto H = hyperbolic_anomaly(mean_anomaly_t{M}, e);

  return rotate_from_orbit_plane(
      cartesian_elements_hyperbola_orbit_plane<num_t, arr_t>(p, e, H, kappa2),
      i, Omega, g);
}

template <numm::number_type num_t, typename arr_t>
constexpr auto
cartesian_elements_parabola_orbit_plane(num_t p, parabolic_sigma_t<num_t> sigma,
                                        num_t kappa2) {
  num_t sigma2 = sigma * sigma,
        p_dot_sigma =
            num_t{2.0} * std::sqrt(kappa2 / p) / (num_t{1.0} + sigma2);

  num_t X = p * (num_t{1.0} - sigma2) / num_t{2.0}, Y = p * sigma,
        dot_X = -sigma * p_dot_sigma, dot_Y = p_dot_sigma;

  return arr_t{X, Y, num_t{0.0}, dot_X, dot_Y, num_t{0.0}};
}

template <numm::number_type num_t, typename arr_t>
constexpr auto
cartesian_elements_parabola(const arr_t &ke,
                            num_t kappa2) // ke: p, 1.0, i, Omega, g, M
{
  num_t p = ke[0], i = ke[2], Omega = ke[3], g = ke[4], M = ke[5];
  auto sigma = parabolic_sigma(mean_anomaly_t{M});

  return rotate_from_orbit_plane(
      cartesian_elements_parabola_orbit_plane<num_t, arr_t>(p, sigma, kappa2),
      i, Omega, g);
}

template <numm::number_type num_t, typename arr_t>
constexpr auto cartesian_elements(const arr_t &ke,
                                  num_t kappa2) // ke: p, e, i, Omega, g, M
{
  if (ke[0] < num_t{0.0} or ke[1] < num_t{0.0})
    return unknown_res<num_t, arr_t>();
  else if (ke[1] >= num_t{0.0} and ke[1] + num_t{8.0} < num_t{9.0})
    return cartesian_elements_ellipse(ke, kappa2);
  else if (ke[1] + num_t{8.0} > num_t{9.0})
    return cartesian_elements_hyperbola(ke, kappa2);
  else
    return cartesian_elements_parabola(ke, kappa2);
}

} // namespace astro

#endif // ASTRO_KEPLER_HPP
