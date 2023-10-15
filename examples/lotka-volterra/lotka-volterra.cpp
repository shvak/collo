#include <chrono>
#include <cmath>
#include <collo/lobatto.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>

constexpr std::size_t so = 2; // system order
constexpr std::size_t ms = 8; // method stages
using sv_t = collo::state_vector_t<double, so>;

struct LotkaVolterra {
  static std::size_t cnt;

  sv_t operator()(double, const auto &y) const {
    ++cnt;
    return {y[0] * (2.0 - y[1]), y[1] * (y[0] - 1.0)};
  }

  static auto count() { return std::exchange(cnt, 0); }

  static double energy(const auto &y) {
    return y[0] - log(y[0]) + y[1] - 2.0 * log(y[1]);
  }
};

std::size_t LotkaVolterra::cnt = 0;

int main(int, char **argv) {
  std::filesystem::current_path(std::filesystem::path{argv[0]}.parent_path());

  std::size_t maxsteps = 500;
  double h = 0.1;

  auto lobatto =
      collo::make_Lobatto<double, so, ms, LotkaVolterra>({1.0, 1.0}, 0.0, h);

  std::size_t iternum = 0;
  double E_start = LotkaVolterra::energy(lobatto.state());

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
                  (LotkaVolterra::energy(lobatto.state()) - E_start) / E_start);
  }
  auto end = std::chrono::steady_clock::now();

  fout.close();

  std::cout << lobatto.steps() << " steps finished\n";
  std::cout << "dE / E = "
            << (LotkaVolterra::energy(lobatto.state()) - E_start) / E_start
            << '\n';
  std::cout << "iterations = " << iternum << '\n';
  std::cout << "RHS invocations = " << LotkaVolterra::count() << '\n';
  std::cout << "time elapsed "
            << std::chrono::duration<double>(end - start).count() << '\n';

  return 0;
}
