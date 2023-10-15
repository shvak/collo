#ifndef LACE_VECTOR_HPP
#define LACE_VECTOR_HPP

#include "numm/base.hpp"
#include <array>

namespace lace {

using numm::number_type;

template <number_type num_t, std::size_t n>
struct vector : public std::array<num_t, n> {
  using base_t = std::array<num_t, n>;

  constexpr num_t &operator[](std::size_t k) { return base_t::at(k); }

  constexpr const num_t &operator[](std::size_t k) const {
    return base_t::at(k);
  }

  static constexpr vector basis(std::size_t k) noexcept {
    vector res{};
    for (auto &x : res)
      x = num_t{0.0};
    res[k] = {1.0};
    return res;
  }

  constexpr base_t to_array() const {
    base_t res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = base_t::at(i);
    return res;
  }

  constexpr vector abs() const {
    vector res{*this};
    for (auto &x : res)
      x = x < num_t{0.0} ? -x : x;
    return res;
  }

  constexpr num_t linorm() const {
    num_t res{0.0};
    for (auto x : this->abs())
      res = std::max(res, x);
    return res;
  }

  constexpr num_t l1norm() const {
    num_t res{0.0};
    for (auto x : this->abs())
      res += x;
    return res;
  }

  constexpr num_t squared_l2norm() const {
    num_t res{0.0};
    for (auto x : *this)
      res += x * x;
    return res;
  }

  constexpr num_t l2norm() const { return numm::sqrt(squared_l2norm()); }

  constexpr num_t norm() const { return l2norm(); }

  constexpr vector operator+(const vector &v) const {
    vector res{*this};
    for (std::size_t i = 0; i < n; ++i)
      res[i] += v[i];
    return res;
  }

  constexpr vector operator+(num_t num) const {
    vector res{*this};
    for (auto &x : res)
      x += num;
    return res;
  }

  constexpr vector operator-(const vector &v) const {
    vector res{*this};
    for (std::size_t i = 0; i < n; ++i)
      res[i] -= v[i];
    return res;
  }

  constexpr vector operator-(num_t num) const {
    vector res{*this};
    for (auto &x : res)
      x -= num;
    return res;
  }

  constexpr vector operator*(num_t num) const {
    vector res{*this};
    for (auto &x : res)
      x *= num;
    return res;
  }

  friend constexpr vector operator*(num_t num, const vector &v) {
    return v * num;
  }

  constexpr vector operator/(num_t num) const {
    vector res{*this};
    for (auto &x : res)
      x /= num;
    return res;
  }

  constexpr num_t operator*(const vector &v) const // dot product
  {
    num_t res{0.0};
    for (std::size_t i = 0; i < n; ++i)
      res += base_t::at(i) * v[i];
    return res;
  }
};

} // namespace lace

#endif // LACE_VECTOR_HPP
