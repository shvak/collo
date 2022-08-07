#ifndef ASTRO_BASE_HPP
#define ASTRO_BASE_HPP

#include <numm/base.hpp>
#include <type_traits>
#include <utility>

namespace astro {

enum class phys_val_tag {
  // semimajor_axis,
  // eccentricity,
  // inclination,
  // longitude_ascending_node,
  // argument_periapsis,
  true_anomaly,
  eccentric_anomaly,
  mean_anomaly,
  d_mean_anomaly,
  hyperbolic_anomaly,
  parabolic_sigma
  // mean_motion
};

template <numm::number_type num_t, phys_val_tag pvt>
struct phys_val /* : public std::pair<num_t, phys_val_tag>*/
{
  // using base_t = std::pair<num_t, phys_val_tag>;
  using val_t = num_t;
  phys_val_tag tag = pvt;
  val_t val;

  phys_val(num_t num) noexcept
      : val{num} // base_t{num, pvt}
  {}

  operator val_t() const noexcept {
    return val; // base_t::first;
  }
};

/*template<numm::number_type num_t>
using semimajor_axis_t = phys_val<num_t, phys_val_tag::semimajor_axis>;

template<numm::number_type num_t>
using eccentricity_t = phys_val<num_t, phys_val_tag::eccentricity>;

template<numm::number_type num_t>
using inclination_t = phys_val<num_t, phys_val_tag::inclination>;

template<numm::number_type num_t>
using longitude_ascending_node_t = phys_val<num_t,
phys_val_tag::longitude_ascending_node>;

template<numm::number_type num_t>
using argument_periapsis_t = phys_val<num_t,
phys_val_tag::argument_periapsis>;*/

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

/*template<numm::number_type num_t>
using mean_motion_t = phys_val<num_t, phys_val_tag::mean_motion>;*/

} // namespace astro

#endif // ASTRO_BASE_HPP
