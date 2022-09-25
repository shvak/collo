#include <astro/kepler.hpp>
#include <chrono>
#include <cmath>
#include <collo/lobatto.hpp>
#include <fmt/ostream.h>
#include <fstream>
#include <iostream>
#include <numbers>

template <typename T>
requires std::is_base_of_v<Eigen::DenseBase<T>, T>
struct fmt::formatter<T> : ostream_formatter {
};

constexpr std::size_t so = 6; // system order
constexpr std::size_t ms = 8; // method stages
using sv_t = collo::state_vector_t<double, so>;
using std::numbers::pi;

struct Kepler {
  static double GM;
  static std::size_t cnt;

  sv_t operator()(double, const auto &y) {
    ++cnt;
    double r2 = y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
    double r = std::sqrt(r2);
    double alpha = -GM / r / r2;
    return {y[3], y[4], y[5], y[0] * alpha, y[1] * alpha, y[2] * alpha};
  }

  static auto count() { return std::exchange(cnt, 0); }

  static double energy(const auto &y) {
    double T = (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]) / 2;
    double r = std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
    double P = -GM / r;
    return T + P;
  }
};

double Kepler::GM = 4 * pi * pi;
std::size_t Kepler::cnt = 0;

int main() {
  std::ifstream fin("start.dat");
  sv_t y0;
  fin >> y0[0] >> y0[1] >> y0[2] >> y0[3] >> y0[4] >>
      y0[5]; // a, e, i, Omega, g, M
  double t0, h;
  fin >> t0;
  fin >> h;
  std::size_t maxsteps;
  fin >> maxsteps;
  fin.close();

  y0[0] *= 1.0 - y0[1] * y0[1]; // from a to p
  y0 = astro::cartesian_elements(y0, Kepler::GM);

  auto lobatto = collo::make_Lobatto<double, so, ms, Kepler>(y0, t0, h);

  std::size_t iternum = 0;
  double E_start = Kepler::energy(lobatto.state());

  std::ofstream fout("out.dat");
  fout << fmt::format("{}\t{}\t{:.17g}\n", /*time*/ 0.0,
                      lobatto.state().transpose(), /*energy*/ 0.0);

  auto start = std::chrono::steady_clock::now();
  while (lobatto.steps() < maxsteps) {
    iternum += lobatto.do_step().iternum();
    fout << fmt::format("{}\t{}\t{:.17g}\n", lobatto.steps() * h,
                        lobatto.state().transpose(),
                        (Kepler::energy(lobatto.state()) - E_start) / E_start);
  }
  auto end = std::chrono::steady_clock::now();

  fout.close();

  std::cout << fmt::format("{} steps finished\n", lobatto.steps());
  std::cout << fmt::format(
      "dE / E = {}\n", (Kepler::energy(lobatto.state()) - E_start) / E_start);
  std::cout << fmt::format("iterations = {}\n", iternum);
  std::cout << fmt::format("RHS invocations = {}\n", Kepler::count());
  std::cout << fmt::format("time elapsed {}",
                           std::chrono::duration<double>(end - start).count());

  return 0;
}
