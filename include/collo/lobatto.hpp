#ifndef COLLO_LOBATTO_HPP
#define COLLO_LOBATTO_HPP

#include "collocation.hpp"

namespace collo {

template <numm::number_type num_t, std::size_t system_order,
          std::size_t method_stage>
struct Lobatto_base
    : protected Collocation_base<num_t, system_order, method_stage> {

protected:
  static constexpr auto time_nodes =
      numm::roots_ilegendre_sh<method_stage - 1, num_t>();
};

template <Pred p, numm::number_type num_t, std::size_t system_order,
          std::size_t method_stage, std::size_t param>
using Pred_Lobatto_base_t = Pred_select<p, num_t, system_order, method_stage,
                                        Lobatto_base, param>::type;

// template <numm::number_type num_t, std::size_t system_order,
//           std::size_t method_stage, typename rhs_t, Pred p = Pred::Poly,
//           std::size_t param = 3,
//           typename sv_t = state_vector_t<num_t, system_order>>
// constexpr auto make_Lobatto(const sv_t &y0, num_t t0, num_t h, rhs_t rhs
//                             /* std::decay_t<rhs_t> rhs */) {
//   return Collocation<
//       num_t, std::decay_t<rhs_t>,
//       Pred_Lobatto_base_t<p, num_t, system_order, method_stage, param>>(y0,
//       t0,
//                                                                         h,
//                                                                         rhs);
// }

// template <numm::number_type num_t, std::size_t system_order,
//           std::size_t method_stage, typename rhs_t, Pred p = Pred::Poly,
//           std::size_t param = 3,
//           typename sv_t = state_vector_t<num_t, system_order>>
// constexpr auto make_Lobatto(const sv_t &y0, num_t t0, num_t h,
//                             rhs_t &&rhs = rhs_t{}) {
//   return Collocation<
//       num_t, rhs_t,
//       Pred_Lobatto_base_t<p, num_t, system_order, method_stage, param>>(
//       y0, t0, h, std::forward<rhs_t>(rhs));
// }

template <numm::number_type num_t, std::size_t system_order,
          std::size_t method_stage, typename rhs_t, Pred p = Pred::Poly,
          std::size_t param = 3,
          typename sv_t = state_vector_t<num_t, system_order>>
constexpr auto make_Lobatto(sv_t &&y0, auto &&t0, auto &&h,
                            rhs_t &&rhs = rhs_t{}) {
  return Collocation<
      num_t, rhs_t,
      Pred_Lobatto_base_t<p, num_t, system_order, method_stage, param>>(
      std::forward<sv_t>(y0), std::forward<num_t>(t0), std::forward<num_t>(h),
      std::forward<rhs_t>(rhs));
}

} // namespace collo

#endif // COLLO_LOBATTO_HPP
