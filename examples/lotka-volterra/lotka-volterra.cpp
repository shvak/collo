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
  fout << fmt::format("{}\t{}\t{:.17g}\n", /*time*/ 0.0,
                      lobatto.state().transpose(), /*energy*/ 0.0);

  auto start = std::chrono::steady_clock::now();
  while (lobatto.steps() < maxsteps) {
    iternum += lobatto.do_step().iternum();
    fout << fmt::format(
        "{}\t{}\t{:.17g}\n", lobatto.steps() * h, lobatto.state().transpose(),
        (LotkaVolterra::energy(lobatto.state()) - E_start) / E_start);
  }
  auto end = std::chrono::steady_clock::now();

  fout.close();

  std::cout << fmt::format("{} steps finished\n", lobatto.steps());
  std::cout << fmt::format("dE / E = {}\n",
                           (LotkaVolterra::energy(lobatto.state()) - E_start) /
                               E_start);
  std::cout << fmt::format("iterations = {}\n", iternum);
  std::cout << fmt::format("RHS invocations = {}\n", LotkaVolterra::count());
  std::cout << fmt::format("time elapsed {}",
                           std::chrono::duration<double>(end - start).count());

  return 0;
}
