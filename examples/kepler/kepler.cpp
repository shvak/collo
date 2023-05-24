#include <astro/kepler.hpp>
#include <chrono>
#include <cmath>
#include <collo/lobatto.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>

constexpr std::size_t so = 6; // system order
constexpr std::size_t ms = 8; // method stages
using sv_t = collo::state_vector_t<double, so>;
using std::numbers::pi;

struct Kepler {
  static constexpr double GM = 4 * pi * pi;
  static std::size_t cnt;

  sv_t operator()(double, const auto &y) const {
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

std::size_t Kepler::cnt = 0;

int main(int, char **argv) {
  std::filesystem::current_path(std::filesystem::path{argv[0]}.parent_path());

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
  fout.precision(15);
  auto write_to_file = [&fout](double time, const sv_t &y, double energy) {
    fout << time << '\t' << y.transpose() << '\t' << energy << '\n';
  };
  write_to_file(/*time*/ 0., lobatto.state(), /*energy*/ 0.);

  auto start = std::chrono::steady_clock::now();
  while (lobatto.steps() < maxsteps) {
    iternum += lobatto.do_step().iternum();
    write_to_file(lobatto.time(), lobatto.state(),
                  (Kepler::energy(lobatto.state()) - E_start) / E_start);
  }
  auto end = std::chrono::steady_clock::now();

  fout.close();

  std::cout << lobatto.steps() << " steps finished\n";
  std::cout << "dE / E = "
            << (Kepler::energy(lobatto.state()) - E_start) / E_start << '\n';
  std::cout << "iterations = " << iternum << '\n';
  std::cout << "RHS invocations = " << Kepler::count() << '\n';
  std::cout << "time elapsed "
            << std::chrono::duration<double>(end - start).count() << '\n';

  return 0;
}
