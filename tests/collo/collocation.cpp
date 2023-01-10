#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "collo/lobatto.hpp"
#include <doctest/doctest.h>

using sv_t = collo::state_vector_t<double, 2>;

struct copy_test_rhs {
  static std::size_t cnt;

  copy_test_rhs() = default;
  copy_test_rhs(const copy_test_rhs &) { ++cnt; }
  copy_test_rhs(copy_test_rhs &&) = default;

  sv_t operator()(double, const auto &y) const { return {y[1], -y[0]}; }
};

std::size_t copy_test_rhs::cnt = 0;

TEST_CASE("collo/collocation.hpp") {
  auto l = collo::make_Lobatto<double, 2, 2, copy_test_rhs>(sv_t{1.0, 0.0}, 0.0,
                                                            1.0);
  while (l.do_step().steps() < 100)
    ;
  CHECK(copy_test_rhs::cnt == 1);
}
