#ifndef LACE_VECTOR_HPP
#define LACE_VECTOR_HPP

#include "numm/base.hpp"
#include <array>

namespace lace {

using numm::number_type;

template <number_type num_t, std::size_t n>
struct vector : public std::array<num_t, n> {
  using base_t = std::array<num_t, n>;

  consteval num_t &operator[](std::size_t k) { return base_t::at(k); }

  consteval const num_t &operator[](std::size_t k) const {
    return base_t::at(k);
  }

  static consteval vector basis(std::size_t k) {
    vector res;
    for (auto &x : res)
      x = num_t{0.0};
    res[k] = {1.0};
    return res;
  }

  consteval base_t to_array() const {
    base_t res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = base_t::at(i);
    return res;
  }

  consteval vector abs() const {
    vector res{*this};
    for (auto &x : res)
      x = x < num_t{0.0} ? -x : x;
    return res;
  }

  consteval num_t linorm() const {
    num_t res{0.0};
    for (auto x : this->abs())
      res = std::max(res, x);
    return res;
  }

  consteval num_t l1norm() const {
    num_t res{0.0};
    for (auto x : this->abs())
      res += x;
    return res;
  }

  consteval num_t squared_l2norm() const {
    num_t res{0.0};
    for (auto x : *this)
      res += x * x;
    return res;
  }

  consteval num_t l2norm() const {
    return numm::sqrt(squared_l2norm());
  }

  consteval num_t norm() const { return l2norm(); }

  consteval vector operator+(const vector &v) const {
    vector res{*this};
    for (std::size_t i = 0; i < n; ++i)
      res[i] += v[i];
    return res;
  }

  consteval vector operator+(num_t num) const {
    vector res{*this};
    for (auto &x : res)
      x += num;
    return res;
  }

  consteval vector operator-(const vector &v) const {
    vector res{*this};
    for (std::size_t i = 0; i < n; ++i)
      res[i] -= v[i];
    return res;
  }

  consteval vector operator-(num_t num) const {
    vector res{*this};
    for (auto &x : res)
      x -= num;
    return res;
  }

  consteval vector operator*(num_t num) const {
    vector res{*this};
    for (auto &x : res)
      x *= num;
    return res;
  }

  friend consteval vector operator*(num_t num, const vector &v) {
    return v * num;
  }

  consteval vector operator/(num_t num) const {
    vector res{*this};
    for (auto &x : res)
      x /= num;
    return res;
  }

  consteval num_t operator*(const vector &v) const // dot product
  {
    num_t res{0.0};
    for (std::size_t i = 0; i < n; ++i)
      res += base_t::at(i) * v[i];
    return res;
  }
};

} // namespace lace

#endif // LACE_VECTOR_HPP
