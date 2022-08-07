#ifndef COLLO_GAUSS_HPP
#define COLLO_GAUSS_HPP

#include "collocation.hpp"

namespace collo {

template <numm::number_type num_t, std::size_t system_order,
          std::size_t method_stage>
struct Gauss_base
    : protected Collocation_base<num_t, system_order, method_stage> {

protected:
  static constexpr auto time_nodes =
      numm::roots_legendre_sh<method_stage, num_t>();
};

template <Pred p, numm::number_type num_t, std::size_t system_order,
          std::size_t method_stage, std::size_t param>
using Pred_Gauss_base_t =
    Pred_select<p, num_t, system_order, method_stage, Gauss_base, param>::type;

template <numm::number_type num_t, std::size_t system_order,
          std::size_t method_stage, typename rhs_t, Pred p = Pred::Poly,
          std::size_t param = 3,
          typename sv_t = state_vector_t<num_t, system_order>>
constexpr auto make_Gauss(const sv_t &y0, num_t t0, num_t h,
                          std::decay_t<rhs_t> rhs) {
  return Collocation<
      num_t, std::decay_t<rhs_t>,
      Pred_Gauss_base_t<p, num_t, system_order, method_stage, param>>(y0, t0, h,
                                                                      rhs);
}

template <numm::number_type num_t, std::size_t system_order,
          std::size_t method_stage, typename rhs_t, Pred p = Pred::Poly,
          std::size_t param = 3,
          typename sv_t = state_vector_t<num_t, system_order>>
constexpr auto make_Gauss(const sv_t &y0, num_t t0, num_t h,
                          rhs_t &&rhs = rhs_t{}) {
  return Collocation<
      num_t, rhs_t,
      Pred_Gauss_base_t<p, num_t, system_order, method_stage, param>>(
      y0, t0, h, std::forward<rhs_t>(rhs));
}

template <numm::number_type num_t, std::size_t system_order,
          std::size_t method_stage, typename rhs_t, Pred p = Pred::Poly,
          std::size_t param = 3,
          typename sv_t = state_vector_t<num_t, system_order>>
constexpr auto make_Gauss(sv_t &&y0, num_t &&t0, num_t &&h,
                          rhs_t &&rhs = rhs_t{}) {
  return Collocation<
      num_t, rhs_t,
      Pred_Gauss_base_t<p, num_t, system_order, method_stage, param>>(
      std::forward<sv_t>(y0), std::forward<num_t>(t0), std::forward<num_t>(h),
      std::forward<rhs_t>(rhs));
}

} // namespace collo

#endif // COLLO_GAUSS_HPP
