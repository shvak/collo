#ifndef LACE_MATRIX_HPP
#define LACE_MATRIX_HPP

#include <lace/vector.hpp>
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

  consteval vec_t col(std::size_t k) const {
    vec_t res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = base_t::at(i)[k];
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

  consteval square_matrix operator*(const num_t &num) const {
    square_matrix res;
    for (std::size_t i = 0; i < n; ++i)
      res[i] = row(i) * num;
    return res;
  }

  friend consteval square_matrix operator*(const num_t &num,
                                           const square_matrix &m) {
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

  consteval auto PLU() const {
    auto swap =
        [](const square_matrix &m, std::size_t i, std::size_t j) consteval {
      square_matrix res{m};
      vec_t v = res.row(i);
      res.row(i) = res.row(j);
      res.row(j) = v;
      return res;
    };
    auto subcol = [](const square_matrix &m, std::size_t i) consteval {
      vec_t res{m.col(i)};
      for (std::size_t k = 0; k < i; ++k)
        res[k] = 0;
      return res;
    };
    auto pivot = [](const vec_t &v) consteval {
      std::size_t res = 0;
      vec_t abs_v = v.abs();
      for (std::size_t i = 0; i < n; ++i)
        if (abs_v[res] < abs_v[i])
          res = i;
      return res;
    };
    auto anti = [](const vec_t &v, std::size_t k) consteval {
      vec_t res = vec_t::basis(k) - v / v[k];
      res[k] = num_t{1.0} / v[k];
      return res;
    };
    auto upd =
        [](const square_matrix &m, const vec_t &v, std::size_t k) consteval {
      square_matrix el{identity()};
      el.row(k) = v;
      return m * el;
    };
    square_matrix E = identity();
    square_matrix P{E}, Lt{E}, Lti{E}, U{*this}, Ui{E};
    for (std::size_t k = 0; k < n; ++k) {
      vec_t l = subcol(U, k);
      std::size_t p = pivot(l);
      if (p != k) {
        P = swap(E, k, p) * P;
        U = swap(U, k, p);
      }
      l = subcol(U, k);
      Lt.row(k) = l;
      vec_t li = anti(l, k);
      Lti = upd(Lti, li, k);
      for (std::size_t j = k + 1; j < n; ++j)
        U.row(j) = U.row(j) + U.row(k) * li[j];
      U.row(k) = U.row(k) / l[k];
      Ui = upd(Ui, anti(U.row(k), k), k);
    }
    return std::make_tuple(P, Lt.transpose(), U, Lti.transpose(), Ui);
  }

  consteval square_matrix inverse() const {
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
