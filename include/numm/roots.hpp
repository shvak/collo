#ifndef NUMM_ROOTS_HPP
#define NUMM_ROOTS_HPP

#include "base.hpp"

namespace numm {

template <typename Func, number_type num_t = double>
constexpr num_t fixed_point(const num_t &start, Func f,
                            const num_t delta = num_t{0.1},
                            const num_t shift_from_zero = num_t{2}) {
  num_t tau = 2 * delta / (f(start - delta) - f(start + delta));
  auto s = [tau, f](num_t x) { return x + tau * f(x); };
  num_t prev{start}, curr{s(start)};
  while (curr - prev + shift_from_zero != shift_from_zero) {
    prev = curr;
    curr = s(curr);
  }
  return curr;
}

template <typename Func, typename Deriv, number_type num_t = double>
constexpr num_t newton(const num_t &start, Func f, const Deriv &df,
                       const num_t shift_from_zero = num_t{2}) {
  num_t prev{start}, curr{start - f(start) / df(start)};
  while (curr - prev + shift_from_zero != shift_from_zero) {
    prev = curr;
    curr = prev - f(prev) / df(prev);
  }
  return curr;
}

} // namespace numm

#endif // NUMM_ROOTS_HPP
