#ifndef COLLO_COLLOCATION_HPP
#define COLLO_COLLOCATION_HPP

#include <Eigen/Dense>
#include <deque>
#include <lace/matrix.hpp>
#include <numm/legendre.hpp>
#include <optional>

namespace collo {

template <numm::number_type num_t, std::size_t system_order>
using state_vector_t = Eigen::Matrix<num_t, system_order, 1>;

template <numm::number_type num_t, std::size_t system_order,
          std::size_t method_stage>
struct Collocation_base {
protected:
  using sv_t = state_vector_t<num_t, system_order>;
  using vector_t = Eigen::Matrix<num_t, method_stage, 1>;
  using matrix_t = Eigen::Matrix<num_t, method_stage, method_stage>;
  using sva_t = Eigen::Matrix<num_t, system_order,
                              method_stage>; // type of state_vector's array
  using lacev = lace::vector<num_t, method_stage>;
  using lacem = lace::square_matrix<num_t, method_stage>;

  constexpr static auto
  shifted_tn(const std::array<num_t, method_stage> &time_nodes,
             num_t shift = 0.0) {
    if (shift == 0)
      return time_nodes;
    auto tn = lacev{time_nodes};
    auto tn_sh = tn + shift;
    return tn_sh.to_array();
  }

  constexpr static auto
  make_nodes_basis(const std::array<num_t, method_stage> &time_nodes, num_t arg,
                   num_t shift) {
    auto nodes_basis = lacem::from_2darray(
        numm::legendre_sh<method_stage, 1>(shifted_tn(time_nodes, shift)));
    auto base_sv_basis = lacev{numm::legendre_sh<method_stage, 1>(arg)};
    auto nodes_basis_mod = nodes_basis - base_sv_basis;
    return nodes_basis_mod.to_2darray();
  }

  constexpr static auto
  make_nodes_basis_left(const std::array<num_t, method_stage> &time_nodes) {
    return make_nodes_basis(time_nodes, 0.0, 0.0);
  }

  constexpr static auto
  make_nodes_basis_right(const std::array<num_t, method_stage> &time_nodes) {
    return make_nodes_basis(time_nodes, 1.0, 0.0);
  }

  constexpr static auto make_pred_nodes_basis_right(
      const std::array<num_t, method_stage> &time_nodes) {
    return make_nodes_basis(time_nodes, 1.0, 1.0);
  }

  constexpr static auto
  make_inv_lsm(const std::array<num_t, method_stage> &time_nodes) {
    auto lsm =
        lacem::from_2darray(numm::dlegendre_sh<method_stage, 1>(time_nodes));
    return lsm.invert().to_1darray();
  }

  static vector_t basis_left(num_t tau) {
    return vector_t{numm::legendre_sh<method_stage, 1>(tau).data()} -
           vector_t{numm::legendre_sh<method_stage, 1>(num_t{1.0}).data()};
  }

  static sv_t result(const sva_t &alphas) {
    using namespace Eigen::indexing;
    return alphas(all, seq(0, last, fix<2>)).eval().rowwise().sum() * 2;
  }

  static num_t distance(const sva_t &first, const sva_t &second) {
    return (first - second).norm();
  }
};

template <numm::number_type num_t, typename rhs_t, typename base_t>
struct Collocation : protected base_t {
private:
  using sv_t = typename base_t::sv_t;
  using sva_t = typename base_t::sva_t;
  using base_t::basis_left;
  using base_t::distance;
  using base_t::inv_lsm;
  using base_t::make_alphas;
  using base_t::node_basis_left;
  using base_t::node_basis_right;
  using base_t::result;
  using base_t::rhs_invoke;
  using base_t::save_alphas;
  using base_t::tp;

  sv_t y;
  num_t t0;
  num_t h;
  std::size_t steps_num{};
  rhs_t rhs;
  std::size_t iter_num{};
  sva_t alphas{};

  sv_t sv_point(std::size_t i) const { return y + alphas * node_basis_left(i); }

public:
  Collocation(auto &&y0, auto &&t0_, auto &&h_, auto &&rhs_)
      : y{std::forward<decltype(y0)>(y0)}, t0{std::forward<decltype(t0_)>(t0_)},
        h{std::forward<decltype(h_)>(h_)},
        rhs{std::forward<decltype(rhs_)>(rhs_)} {}

  auto time() const { return t0 + steps_num * h; }

  auto time_point(std::size_t i) const { return tp(time(), h, i); }

  const auto &state() const { return y; }

  auto iternum() const { return iter_num; }

  auto steps() const { return steps_num; }

  const rhs_t &force() const { return rhs; }

  rhs_t &force() { return rhs; }

  sv_t poly(num_t tau) const { return y + alphas * basis_left(tau); }

  sv_t poly_node(std::size_t i) const {
    return y + alphas * node_basis_right(i);
  }

  Collocation &do_step() {
    iter_num = 0;
    num_t dist;

    alphas = make_alphas(alphas, y, time(), h, rhs);

    do {
      auto alphas_prev = alphas;
      sva_t f;
      for (std::size_t i = 0; i < sva_t::ColsAtCompileTime; ++i)
        f.col(i) =
            rhs_invoke(rhs, time_point(i), std::make_optional(i), sv_point(i));
      alphas = f * inv_lsm() * h;
      dist = distance(alphas, alphas_prev);
      ++iter_num;
    } while (dist + num_t{2.0} != num_t{2.0} && iter_num < 100);

    save_alphas(alphas);

    y += result(alphas);
    ++steps_num;
    return *this;
  }
};

template <typename base_t> struct Predictor_Base : protected base_t {
private:
  using matrix_t = typename base_t::matrix_t;
  using vector_t = typename base_t::vector_t;

protected:
  static const auto &inv_lsm() {
    // auto transpose: array is row major, matrix is col major
    static matrix_t inv_lsm =
        matrix_t{base_t::make_inv_lsm(base_t::time_nodes).data()};
    return inv_lsm;
  }

  static auto node_basis_left(std::size_t i) {
    static constexpr auto nodes_basis =
        base_t::make_nodes_basis_left(base_t::time_nodes);
    return Eigen::Map<const vector_t>(nodes_basis[i].data());
  }

  static auto node_basis_right(std::size_t i) {
    constexpr static auto nodes_basis =
        base_t::make_nodes_basis_right(base_t::time_nodes);
    return Eigen::Map<const vector_t>(nodes_basis[i].data());
  }

  constexpr static auto time_node(std::size_t i) {
    return base_t::time_nodes[i];
  }

  static auto tp(auto t, auto h, std::size_t i) { return t + time_node(i) * h; }

  static auto dt(auto h, std::size_t i) {
    if (i == 0)
      return time_node(0) * h;
    else
      return (time_node(i) - time_node(i - 1)) * h;
  }

  static auto rhs_invoke(const auto &rhs, auto t, std::optional<std::size_t> i,
                         const auto &y) {
    if constexpr (std::is_invocable_v<decltype(rhs), decltype(t),
                                      std::optional<std::size_t>,
                                      typename base_t::sv_t>)
      return rhs(t, i, y);
    else
      return rhs(t, y);
  }
};

template <typename base_t>
struct Predictor_Zero : protected Predictor_Base<base_t> {
protected:
  using sv_t = typename base_t::sv_t;
  using sva_t = typename base_t::sva_t;

  auto make_alphas(const sva_t &, const sv_t &, auto, auto,
                   const auto &) const {
    return sva_t::Zero();
  }

  void save_alphas(const sva_t &) {}
};

template <typename base_t>
struct Predictor_Simple : protected Predictor_Base<base_t> {
protected:
  using sv_t = typename base_t::sv_t;
  using sva_t = typename base_t::sva_t;
  using Predictor_Base<base_t>::inv_lsm;
  using Predictor_Base<base_t>::rhs_invoke;

  sva_t make_alphas(const sva_t &, const sv_t &y, auto t, auto h,
                    const auto &rhs) const {
    sva_t f;
    for (std::size_t i = 0; i < sva_t::ColsAtCompileTime; ++i)
      f.col(i) = rhs_invoke(rhs, t, std::nullopt, y);
    return f * inv_lsm() * h;
  }

  void save_alphas(const sva_t &) {}
};

template <typename base_t>
struct Predictor_PrevStep : protected Predictor_Base<base_t> {
protected:
  using sv_t = typename base_t::sv_t;
  using sva_t = typename base_t::sva_t;

  auto make_alphas(const sva_t &alphas, const sv_t &, auto, auto,
                   const auto &) const {
    return alphas;
  }

  void save_alphas(const sva_t &) {}
};

template <typename base_t>
struct Predictor_Euler : protected Predictor_Base<base_t> {
protected:
  using sv_t = typename base_t::sv_t;
  using sva_t = typename base_t::sva_t;
  using Predictor_Base<base_t>::tp;
  using Predictor_Base<base_t>::dt;
  using Predictor_Base<base_t>::inv_lsm;
  using Predictor_Base<base_t>::rhs_invoke;

private:
  static sv_t euler(const sv_t &y, auto t, std::optional<std::size_t> i, auto h,
                    const auto &rhs) {
    if (h == 0)
      return y;
    return y + rhs_invoke(rhs, t, i, y) * h;
  }

protected:
  sva_t make_alphas(const sva_t &, const sv_t &y, auto t, auto h,
                    const auto &rhs) const {
    sva_t f;
    auto yt = euler(y, t, std::nullopt, dt(h, 0), rhs);
    f.col(0) =
        rhs_invoke(rhs, tp(t, h, 0), std::make_optional<std::size_t>(0), yt);
    for (std::size_t i = 1; i < sva_t::ColsAtCompileTime; ++i) {
      yt = euler(yt, tp(t, h, i - 1), std::make_optional(i - 1), dt(h, i), rhs);
      f.col(i) = rhs_invoke(rhs, tp(t, h, i), std::make_optional(i), yt);
    }
    return f * inv_lsm() * h;
  }

  void save_alphas(const sva_t &) {}
};

template <typename base_t>
struct Predictor_RK4 : protected Predictor_Base<base_t> {
protected:
  using sv_t = typename base_t::sv_t;
  using sva_t = typename base_t::sva_t;
  using Predictor_Base<base_t>::tp;
  using Predictor_Base<base_t>::dt;
  using Predictor_Base<base_t>::inv_lsm;
  using Predictor_Base<base_t>::rhs_invoke;

private:
  static sv_t rk4(const sv_t &y, auto t, std::optional<std::size_t> i, auto h,
                  const auto &rhs) {
    if (h == 0)
      return y;
    auto k1 = rhs_invoke(rhs, t, i, y);
    auto k2 = rhs_invoke(rhs, t + h / 2, i, y + k1 / 2);
    auto k3 = rhs_invoke(rhs, t + h / 2, i, y + k2 / 2);
    if (i)
      ++(i.value());
    auto k4 = rhs_invoke(rhs, t + h, i, y + k3);
    return y + (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6) * h;
  }

protected:
  sva_t make_alphas(const sva_t &, const sv_t &y, auto t, auto h,
                    const auto &rhs) const {
    sva_t f;
    auto yt = rk4(y, t, std::nullopt, dt(h, 0), rhs);
    f.col(0) =
        rhs_invoke(rhs, tp(t, h, 0), std::make_optional<std::size_t>(0), yt);
    for (std::size_t i = 1; i < sva_t::ColsAtCompileTime; ++i) {
      yt = rk4(yt, tp(t, h, i - 1), std::make_optional(i - 1), dt(h, i), rhs);
      f.col(i) = rhs_invoke(rhs, tp(t, h, i), std::make_optional(i), yt);
    }
    return f * inv_lsm() * h;
  }

  void save_alphas(const sva_t &) {}
};

template <typename base_t>
struct Predictor_Poly : protected Predictor_Base<base_t> {
protected:
  using sv_t = typename base_t::sv_t;
  using sva_t = typename base_t::sva_t;
  using vector_t = typename base_t::vector_t;
  using Predictor_Base<base_t>::tp;
  using Predictor_Base<base_t>::inv_lsm;
  using Predictor_Base<base_t>::rhs_invoke;

private:
  static auto pred_node_basis_right(std::size_t i) {
    constexpr static auto nodes_basis =
        base_t::make_pred_nodes_basis_right(base_t::time_nodes);
    return Eigen::Map<const vector_t>(nodes_basis[i].data());
  }

protected:
  sva_t make_alphas(const sva_t &alphas, const sv_t &y, auto t, auto h,
                    const auto &rhs) const {
    sva_t f;
    for (std::size_t i = 0; i < sva_t::ColsAtCompileTime; ++i) {
      sv_t yt = y + alphas * pred_node_basis_right(i);
      f.col(i) = rhs_invoke(rhs, tp(t, h, i), std::make_optional(i), yt);
    }
    return f * inv_lsm() * h;
  }

  void save_alphas(const sva_t &) {}
};

template <typename base_t, std::size_t k>
struct Predictor_BackDiff : protected Predictor_Base<base_t> {
protected:
  using sv_t = typename base_t::sv_t;
  using sva_t = typename base_t::sva_t;

private:
  std::array<std::deque<sva_t>, k> bdiff{};

protected:
  using Predictor_Base<base_t>::tp;
  using Predictor_Base<base_t>::inv_lsm;
  using Predictor_Base<base_t>::node_basis_left;
  using Predictor_Base<base_t>::rhs_invoke;

  // Predictor_BackDiff() :  {}

  sva_t make_alphas(const sva_t &, const sv_t &y, auto t, auto h,
                    const auto &rhs) const {
    auto lim = bdiff.front().size();
    if (lim == 0)
      return sva_t::Zero().eval();
    auto zs = bdiff.front().front();
    for (std::size_t i = 1; i < lim; ++i)
      zs += bdiff[i].front();
    sva_t f;
    for (std::size_t i = 0; i < sva_t::ColsAtCompileTime; ++i) {
      sv_t yt = y + zs.col(i);
      f.col(i) = rhs_invoke(rhs, tp(t, h, i), std::make_optional(i), yt);
    }
    return f * inv_lsm() * h;
  }

  void save_alphas(const sva_t &alphas) {
    sva_t z;
    for (std::size_t i = 0; i < sva_t::ColsAtCompileTime; ++i)
      z.col(i) = alphas * node_basis_left(i);
    bdiff.front().push_front(z);
    std::size_t lim = bdiff.front().size();
    for (std::size_t i = 1; i < std::min(lim, k); ++i)
      bdiff[i].push_front(bdiff[i - 1][0] - bdiff[i - 1][1]);
    if (lim > k)
      for (auto &bd : bdiff)
        bd.pop_back();
  }
};

enum class Pred { Zero, Simple, PrevStep, Euler, RK4, Poly, BackDiff };

// template <Pred p, numm::number_type num_t, std::size_t system_order,
//           std::size_t method_stage,
//           template <numm::number_type, std::size_t, std::size_t>
//           typename base_t,
//           std::size_t param>
// struct Pred_select {
//   using type = std::conditional_t<
//       p == Pred::Zero,
//       Predictor_Zero<base_t<num_t, system_order, method_stage>>,
//       std::conditional_t<
//           p == Pred::Simple,
//           Predictor_Simple<base_t<num_t, system_order, method_stage>>,
//           std::conditional_t<
//               p == Pred::PrevStep,
//               Predictor_PrevStep<base_t<num_t, system_order, method_stage>>,
//               std::conditional_t<
//                   p == Pred::Euler,
//                   Predictor_Euler<base_t<num_t, system_order, method_stage>>,
//                   std::conditional_t<
//                       p == Pred::RK4,
//                       Predictor_RK4<base_t<num_t, system_order,
//                       method_stage>>, std::conditional_t<
//                           p == Pred::Poly,
//                           Predictor_Poly<
//                               base_t<num_t, system_order, method_stage>>,
//                           std::conditional_t<
//                               p == Pred::BackDiff,
//                               Predictor_BackDiff<
//                                   base_t<num_t, system_order, method_stage>,
//                                   param>,
//                               void>>>>>>>;
// };
//
template <Pred p, typename base_t, std::size_t param> struct Pred_select {
  template <Pred V1, Pred V2, typename T>
  struct case_t : std::bool_constant<V1 == V2> {
    using type = T;
  };

  template <typename T> struct default_type : std::true_type {
    using type = T;
  };

  using type = typename std::disjunction<
      case_t<p, Pred::Zero, Predictor_Zero<base_t>>,
      case_t<p, Pred::Simple, Predictor_Simple<base_t>>,
      case_t<p, Pred::PrevStep, Predictor_PrevStep<base_t>>,
      case_t<p, Pred::Euler, Predictor_Euler<base_t>>,
      case_t<p, Pred::RK4, Predictor_RK4<base_t>>,
      case_t<p, Pred::Poly, Predictor_Poly<base_t>>,
      case_t<p, Pred::BackDiff, Predictor_BackDiff<base_t, param>>,
      default_type<void>>::type;
  //
  // using type = std::conditional_t<
  //     p == Pred::Zero, Predictor_Zero<base_t>,
  //     std::conditional_t<
  //         p == Pred::Simple, Predictor_Simple<base_t>,
  //         std::conditional_t<
  //             p == Pred::PrevStep, Predictor_PrevStep<base_t>,
  //             std::conditional_t<
  //                 p == Pred::Euler, Predictor_Euler<base_t>,
  //                 std::conditional_t<
  //                     p == Pred::RK4, Predictor_RK4<base_t>,
  //                     std::conditional_t<
  //                         p == Pred::Poly, Predictor_Poly<base_t>,
  //                         std::conditional_t<p == Pred::BackDiff,
  //                                            Predictor_BackDiff<base_t,
  //                                            param>, void>>>>>>>;
};

} // namespace collo

#endif // COLLO_COLLOCATION_HPP
