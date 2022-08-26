#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "collo/gauss.hpp"
#include <doctest/doctest.h>

using sv_t = collo::state_vector_t<double, 2>;

sv_t f(double, const sv_t &y) { return sv_t{y[1], -y[0]}; }

TEST_CASE("collo/gauss.hpp") {
  // auto f = [](double, const sv_t &y) { return sv_t{y[1], -y[0]}; };
  auto l = collo::make_Gauss<double, 2, 2>(sv_t{1.0, 0.0}, 0.0, 1.0, f);
  CHECK(l.state() == sv_t{1.0, 0.0});
  CHECK(l.steps() == 0);
  CHECK(l.time() == 0.0);
  CHECK(l.iternum() == 0);
  while (l.do_step().steps() < 1)
    ;
  CHECK(l.steps() == 1);
  CHECK(l.time() == 1.0);
}
