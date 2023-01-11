#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "collo/lobatto.hpp"
#include <doctest/doctest.h>

using sv_t = collo::state_vector_t<double, 2>;

struct test_rhs {
  sv_t operator()(double, const auto &y) const { return {y[1], -y[0]}; }
};

TEST_CASE("collo/collocation.hpp: different predictors") {
  auto l_zero = collo::make_Lobatto<double, 2, 2, test_rhs, collo::Pred::Zero>(
      sv_t{1.0, 0.0}, 0.0, 1.0);
  while (l_zero.do_step().steps() < 10)
    ;
  CHECK(l_zero.steps() == 10);
  auto l_simple =
      collo::make_Lobatto<double, 2, 2, test_rhs, collo::Pred::Simple>(
          sv_t{1.0, 0.0}, 0.0, 1.0);
  while (l_simple.do_step().steps() < 10)
    ;
  CHECK(l_simple.steps() == 10);
  auto l_prevstep =
      collo::make_Lobatto<double, 2, 2, test_rhs, collo::Pred::PrevStep>(
          sv_t{1.0, 0.0}, 0.0, 1.0);
  while (l_prevstep.do_step().steps() < 10)
    ;
  CHECK(l_prevstep.steps() == 10);
  auto l_euler =
      collo::make_Lobatto<double, 2, 2, test_rhs, collo::Pred::Euler>(
          sv_t{1.0, 0.0}, 0.0, 1.0);
  while (l_euler.do_step().steps() < 10)
    ;
  CHECK(l_euler.steps() == 10);
  auto l_rk4 = collo::make_Lobatto<double, 2, 2, test_rhs, collo::Pred::RK4>(
      sv_t{1.0, 0.0}, 0.0, 1.0);
  while (l_rk4.do_step().steps() < 10)
    ;
  CHECK(l_rk4.steps() == 10);
  auto l_poly = collo::make_Lobatto<double, 2, 2, test_rhs, collo::Pred::Poly>(
      sv_t{1.0, 0.0}, 0.0, 1.0);
  while (l_poly.do_step().steps() < 10)
    ;
  CHECK(l_poly.steps() == 10);
  auto l_backdiff =
      collo::make_Lobatto<double, 2, 2, test_rhs, collo::Pred::BackDiff>(
          sv_t{1.0, 0.0}, 0.0, 1.0);
  while (l_backdiff.do_step().steps() < 10)
    ;
  CHECK(l_backdiff.steps() == 10);
}

struct copy_test_rhs {
  static std::size_t cnt;

  copy_test_rhs() = default;
  copy_test_rhs(const copy_test_rhs &) { ++cnt; }
  copy_test_rhs(copy_test_rhs &&) = default;

  sv_t operator()(double, const auto &y) const { return {y[1], -y[0]}; }
};

std::size_t copy_test_rhs::cnt = 0;

TEST_CASE("collo/collocation.hpp: copy test for rhs") {
  auto l = collo::make_Lobatto<double, 2, 2, copy_test_rhs>(sv_t{1.0, 0.0}, 0.0,
                                                            1.0);
  while (l.do_step().steps() < 100)
    ;
  CHECK(copy_test_rhs::cnt == 1);
}

struct alt_rhs {
  static std::size_t cnt;

  sv_t operator()(double, std::size_t, const auto &y) const {
    ++cnt;
    return {y[1], -y[0]};
  }
};

std::size_t alt_rhs::cnt = 0;

TEST_CASE("collo/collocation.hpp: alternative rhs form") {
  auto l = collo::make_Lobatto<double, 2, 2, alt_rhs>(sv_t{1.0, 0.0}, 0.0, 1.0);
  while (l.do_step().steps() < 100)
    ;
  CHECK(alt_rhs::cnt > 0);
}
