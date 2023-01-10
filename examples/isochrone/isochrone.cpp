#include <chrono>
#include <cmath>
#include <collo/lobatto.hpp>
#include <filesystem>
#include <fmt/ostream.h>
#include <fstream>
#include <iostream>

template <typename T>
requires std::is_base_of_v<Eigen::DenseBase<T>, T>
struct fmt::formatter<T> : ostream_formatter {
};

constexpr std::size_t so = 4; // system order
constexpr std::size_t ms = 8; // method stages
using sv_t = collo::state_vector_t<double, so>;

struct Isochrone {
  static double GM;
  static double b;
  static double b2;
  static std::size_t cnt;

  sv_t operator()(double, const auto &y) const {
    ++cnt;
    double r2 = y[0] * y[0] + y[1] * y[1];
    double r_mod = std::sqrt(r2 + b2);
    double rb2 = r_mod + b;
    rb2 *= rb2;
    double alpha = GM / r_mod / rb2;
    return {y[2], y[3], -y[0] * alpha, -y[1] * alpha};
  }

  static auto count() { return std::exchange(cnt, 0); }

  static double energy(const auto &y) {
    double T = (y[2] * y[2] + y[3] * y[3]) / 2;
    double r2 = y[0] * y[0] + y[1] * y[1];
    double P = -GM / (b + std::sqrt(r2 + b2));
    return T + P;
  }
};

double Isochrone::GM = 4.49829709181306e-12 * 7.8e9;
double Isochrone::b = 0.15;
double Isochrone::b2 = b * b;
std::size_t Isochrone::cnt = 0;

int main(int, char **argv) {
  std::filesystem::current_path(std::filesystem::path{argv[0]}.parent_path());

  constexpr double vu =
      977.8140491848493; // velocity unit = 1 kpc / 1 My in (km / s)

  std::ifstream fin("start.dat");
  sv_t y0;
  fin >> y0[0] >> y0[1] >> y0[2] >>
      y0[3]; // x in kpc, y in kpc, dot x in km/s, dot y in km/s
  y0[2] /= vu;
  y0[3] /= vu;
  double t0, h;
  fin >> t0;
  fin >> h;
  std::size_t maxsteps;
  fin >> maxsteps;
  fin.close();

  auto lobatto = collo::make_Lobatto<double, so, ms, Isochrone>(y0, t0, h);

  std::size_t iternum = 0;
  double E_start = Isochrone::energy(lobatto.state());

  std::ofstream fout("out.dat");
  fout << fmt::format("{}\t{}\t{:.17g}\n", /*time*/ 0.0,
                      lobatto.state().transpose(), /*energy*/ 0.0);

  auto start = std::chrono::steady_clock::now();
  while (lobatto.steps() < maxsteps) {
    iternum += lobatto.do_step().iternum();
    fout << fmt::format(
        "{}\t{}\t{:.17g}\n", lobatto.steps() * h, lobatto.state().transpose(),
        (Isochrone::energy(lobatto.state()) - E_start) / E_start);
  }
  auto end = std::chrono::steady_clock::now();

  fout.close();

  std::cout << fmt::format("{} steps finished\n", lobatto.steps());
  std::cout << fmt::format("dE / E = {}\n",
                           (Isochrone::energy(lobatto.state()) - E_start) /
                               E_start);
  std::cout << fmt::format("iterations = {}\n", iternum);
  std::cout << fmt::format("RHS invocations = {}\n", Isochrone::count());
  std::cout << fmt::format("time elapsed {}",
                           std::chrono::duration<double>(end - start).count());

  return 0;
}
