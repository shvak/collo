#ifndef LACE_MATRIX_HPP
#define LACE_MATRIX_HPP

#include "lace/vector.hpp"
#include <tuple>

namespace lace {

template <number_type num_t, std::size_t n>
struct square_matrix : std::array<vector<num_t, n>, n> {
  using vec_t = lace::vector<num_t, n>;
  using base_t = std::array<vec_t, n>;

  consteval vec_t &operator[](std::size_t k) { return base_t::at(k); }

  consteval const vec_t &operator[](std::size_t k) const {
    return base_t::at(k);
  }

  static consteval square_matrix identity() {
    square_matrix res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = vec_t::basis(i);
    return res;
  }

  static consteval square_matrix
  from_2darray(const std::array<std::array<num_t, n>, n> &arr) {
    square_matrix res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = vec_t{arr[i]};
    return res;
  }

  consteval auto to_2darray() const {
    std::array<std::array<num_t, n>, n> res;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t o = 0; o < n; ++o)
        res[i][o] = row(i)[o];
    return res;
  }

  consteval auto to_1darray() const {
    std::array<num_t, n * n> res;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t o = 0; o < n; ++o)
        res[i * n + o] = row(i)[o];
    return res;
  }

  consteval vec_t &row(std::size_t k) { return base_t::at(k); }

  consteval const vec_t &row(std::size_t k) const { return base_t::at(k); }

  consteval const vec_t col(std::size_t k) const {
    vec_t res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = row(i)[k];
    return res;
  }

  consteval square_matrix transpose() const {
    square_matrix res;
    for (std::size_t i = 0; i < n; ++i)
      res.row(i) = col(i);
    return res;
  }

  consteval num_t norm_rows() const {
    vec_t norms;
    for (std::size_t i = 0; i < n; ++i)
      norms[i] = row(i).l1norm();
    return norms.linorm();
  }

  consteval num_t norm_cols() const {
    vec_t norms;
    for (std::size_t i = 0; i < n; ++i)
      norms[i] = col(i).l1norm();
    return norms.linorm();
  }

  consteval square_matrix operator+(const vec_t &v) const {
    square_matrix res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = row(i) + v;
    return res;
  }

  consteval square_matrix operator+(const square_matrix &m) const {
    square_matrix res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = row(i) + m[i];
    return res;
  }

  consteval square_matrix operator-(const vec_t &v) const {
    square_matrix res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = row(i) - v;
    return res;
  }

  consteval square_matrix operator-(const square_matrix &m) const {
    square_matrix res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = row(i) - m[i];
    return res;
  }

  consteval square_matrix operator*(num_t num) const {
    square_matrix res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = row(i) * num;
    return res;
  }

  friend consteval square_matrix operator*(num_t num, const square_matrix &m) {
    return m * num;
  }

  consteval vec_t operator*(const vec_t &v) const {
    vec_t res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = row(i) * v;
    return res;
  }

  friend consteval vec_t operator*(const vec_t &v, const square_matrix &m) {
    vec_t res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = v * m.col(i);
    return res;
  }

  consteval square_matrix operator*(const square_matrix &m) const {
    square_matrix res;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        res[i][j] = row(i) * m.col(j);
    return res;
  }

  consteval square_matrix swap(std::size_t i, std::size_t j) const {
    if (i == j)
      return *this;
    square_matrix res{*this};
    auto v = res.row(i);
    res.row(i) = res.row(j);
    res.row(j) = v;
    return res;
  }

  consteval std::size_t find_pivot(std::size_t k) const {
    std::size_t p = k;
    auto v = (*this).col(k).abs();
    for (std::size_t i = k + 1; i < n; ++i)
      if (v[i] > v[p])
        p = i;
    return p;
  }

  consteval auto PLU() const {
    square_matrix L{identity()}, U{*this}, Li{identity()}, Uti{identity()},
        P{identity()};
    for (std::size_t k = 0; k < n; ++k) {
      auto p = U.find_pivot(k);
      P = P.swap(k, p);
      U = U.swap(k, p);
      Li = Li.swap(k, p);
      Li[k][p] = Li[p][k] = num_t{0.0};
      Li[k][k] = Li[p][p] = num_t{1.0};
      L = L.swap(k, p);
      L[k][p] = L[p][k] = num_t{0.0};
      L[k][k] = L[p][p] = num_t{1.0};
      auto d = U[k][k];
      L[k][k] = d;
      auto li = identity();
      li[k][k] = num_t{1.0} / d;
      U.row(k) = U.row(k) / d;
      U[k][k] = num_t{1.0};
      for (std::size_t i = k + 1; i < n; ++i) {
        L[i][k] = U[i][k];
        li[i][k] = -U[i][k] / d;
        U.row(i) = U.row(i) - U.row(k) * U[i][k];
        U[i][k] = num_t{0.0};
      }
      for (std::size_t i = k; i > 0; --i)
        Uti[k][i-1] = -(Uti.row(k) * U.row(i-1));
      Li = li * Li;
    }
    return std::make_tuple(P, L, U, Li, Uti.transpose());
  }

  consteval square_matrix invert() const {
    /*square_matrix im = transpose() * ( num_t{1.0} / norm_rows() / norm_cols()
    ); square_matrix E = identity(); square_matrix err = E - im * (*this);
    while( err.norm_rows() + num_t{2.0} != num_t{2.0} or err.norm_cols() +
    num_t{2.0} != num_t{2.0} )
    {
        im = (E + err + err * err) * im;
        err = E - im * (*this);
    }
    return im;*/
    auto plu = PLU();
    return std::get<4>(plu) * std::get<3>(plu) * std::get<0>(plu);
  }
};

} // namespace lace

#endif // LACE_MATRIX_HPP
