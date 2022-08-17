#ifndef NUMM_ROOTS_HPP
#define NUMM_ROOTS_HPP

#include "base.hpp"

namespace numm {

template <typename Func, number_type num_t = double>
constexpr num_t fixed_point(const num_t &start, Func f,
                            const num_t delta = num_t{0.1},
                            const num_t shift_from_zero = num_t{1.0}) {
  num_t tau = 2 * delta / (f(start + delta) - f(start - delta));
  num_t prev{start}, curr{start - tau * f(start)}, fval{f(curr)};
  while (curr - prev + shift_from_zero != shift_from_zero and
         fval + shift_from_zero != shift_from_zero) {
    prev = curr;
    curr = prev - tau * fval;
    fval = f(curr);
  }
  return curr;
}

template <typename Func, typename Deriv, number_type num_t = double>
constexpr num_t newton(const num_t &start, Func f, Deriv df,
                       const num_t shift_from_zero = num_t{1.0}) {
  num_t prev{start}, curr{start - f(start) / df(start)}, fval{curr};
  while (curr - prev + shift_from_zero != shift_from_zero and
         fval + shift_from_zero != shift_from_zero) {
    prev = curr;
    curr = prev - fval / df(prev);
    fval = f(curr);
  }
  return curr;
}

} // namespace numm

#endif // NUMM_ROOTS_HPP
