#ifndef ASTRO_BASE_HPP
#define ASTRO_BASE_HPP

#include <chrono>
#include <numm/base.hpp>
#include <type_traits>
#include <utility>

namespace astro {

enum class phys_val_tag {
  true_anomaly,
  eccentric_anomaly,
  mean_anomaly,
  d_mean_anomaly,
  hyperbolic_anomaly,
  parabolic_sigma
};

template <numm::number_type num_t, phys_val_tag pvt> struct phys_val {
  using val_t = num_t;
  phys_val_tag tag = pvt;
  val_t val;

  phys_val(num_t num) noexcept : val{num} {}

  operator val_t() const noexcept { return val; }
};

template <numm::number_type num_t>
using true_anomaly_t = phys_val<num_t, phys_val_tag::true_anomaly>;

template <numm::number_type num_t>
using eccentric_anomaly_t = phys_val<num_t, phys_val_tag::eccentric_anomaly>;

template <numm::number_type num_t>
using mean_anomaly_t = phys_val<num_t, phys_val_tag::mean_anomaly>;

template <numm::number_type num_t>
using d_mean_anomaly_t = phys_val<num_t, phys_val_tag::d_mean_anomaly>;

template <numm::number_type num_t>
using hyperbolic_anomaly_t = phys_val<num_t, phys_val_tag::hyperbolic_anomaly>;

template <numm::number_type num_t>
using parabolic_sigma_t = phys_val<num_t, phys_val_tag::parabolic_sigma>;

template <numm::number_type num_t>
constexpr num_t julian_date(
    const std::chrono::year_month_day &date) noexcept {
  int y = static_cast<int>(date.year()),
      m = static_cast<unsigned>(date.month()),
      d = static_cast<unsigned>(date.day());
  int a = (m - 14) / 12;
  num_t jd = (1461 * (y + 4800 + a)) / 4;
  jd += (367 * (m - 2 - 12 * a)) / 12;
  jd -= (3 * ((y + 4900 + a) / 100)) / 4;
  jd += d - 32075 - 0.5;
  return jd;
}

template <numm::number_type num_t>
constexpr num_t julian_date(
    auto &&date,
    const std::chrono::duration<num_t, std::ratio<86400>> &time) noexcept {
  return julian_date<num_t>(std::forward<decltype(date)>(date)) + time.count();
}

template <numm::number_type num_t>
constexpr std::chrono::year_month_day gregorian_from_jd(num_t jd) noexcept {
  int j = jd + .5;
  int f = j + 1401 + (((4 * j + 274277) / 146097) * 3) / 4 - 38;
  int e = 4 * f + 3;
  int g = (e % 1461) / 4;
  int h = 5 * g + 2;
  int d = (h % 153) / 5 + 1;
  int m = (h / 153 + 2) % 12 + 1;
  int y = e / 1461 - 4716 + (14 - m) / 12;
  return std::chrono::day{static_cast<unsigned>(d)} / m / y;
}

template <numm::number_type num_t>
constexpr auto time_from_jd(num_t jd) noexcept {
  double sjd = jd + .5;
  int j = sjd;
  return std::chrono::duration<num_t, std::ratio<86400>>(sjd - j);
}

} // namespace astro

#endif // ASTRO_BASE_HPP
