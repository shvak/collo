#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "collo/lobatto.hpp"
#include <doctest/doctest.h>

TEST_CASE("collo/lobatto.hpp") {
  using sv_t = collo::state_vector_t<double, 2>;
  auto f = [](double, const sv_t &y) { return sv_t{y[1], -y[0]}; };
  auto l = collo::make_Lobatto<double, 2, 2>(sv_t{1.0, 0.0}, 0.0, 1.0, f);
  CHECK(l.state() == sv_t{1.0, 0.0});
  CHECK(l.steps() == 0);
  CHECK(l.time() == 0.0);
  CHECK(l.iternum() == 0);
  while (l.do_step().steps() < 1)
    ;
  CHECK(l.steps() == 1);
  CHECK(l.time() == 1.0);
}
