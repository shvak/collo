# collo

This is the header-only library for collocation integrator for ODEs (ordinary differential equations).

The integrator offers an adjustable order, single-step, implicit method for solving ODEs : $y' = f(t, y)$.
Here $y$ is a n-dimensional vector that is represented by `Eigen::Matrix`.

You can find details about collocation integrators in 
Ernst Hairer, Gerhard Wanner and Christian Lubich “Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations”, 2nd ed. Springer Berlin, Heidelberg, 2006.

## Getting started

You can define $f(t, y)$ as a function object with two arguments: time $t$ and vector $y$.
Use `auto` as a type of second argument or use the helper alias `collo::state_vector_t<some_float_type, n>`.

Include `<collo/lobatto.hpp>` and then you can create an integrator with the function `make_Lobatto`. Example:

    auto lobatto = make_Lobatto<some_float_type, n, s, type_of_f>(y_0, t_0, h);

or

    auto lobatto = make_Lobatto<some_float_type, n, s>(y_0, t_0, h, f);

The `some_float_type` is currently `float`, `double` or `long double`.
Parameter `n` is the dimesion of the vector $y$.
Parameter `s` determines how many intermediate points are used for the calculation of result on every step. This parameter impacts accuracy of integrator $O(h^{2s-2})$.
Arguments `y_0` and `t_0` are initial vector and time. Argument `h` is a time step.

### Methods of the integrator class.

Use `do_step()` to integrate a single step with a fixed step size.

* `steps()` returns the number of calculated steps.
* `time()` returns the current time `t = t_0 + h · steps()`.
* `iternum()` returns the number of iterations of the main implicit method loop on the last step.
* `state()` returns the current state vector $y$.
* `force()` can be used to get access to $f(t, y)$ and change it if needed.
* `poly(some_float_type t)` returns the intermediate interpolated vector $y$ with the argument `t` being between 0 and 1. Other values of argument will give extrapolation with collocation polynomial.
* `poly_node(size_t i)` returns the vector $y$ at collocation node with number `i` ${} \at [0, s)$.

## Compilation

This project is dependent on the library **Eigen** v3.4+ (via `find_package()` in **CMake**). You can get it from [official page](https://eigen.tuxfamily.org/index.php?title=Main_Page).

    git clone https://github.com/shvak/collo.git
    cd collo
    mkdir build
    CXXFLAGS="-O3 -march=native" cmake -DBUILD_TESTING=OFF -S . -B build
    cmake --build build

Note: The integrator performance improves when using the compilation option `-march=native`. This feature is provided by the **Eigen** library and modern CPU vectorization.

Additionally, if you want to use built-in tests, you will need the library **doctest** (via `find_package()` in **CMake**). Then:

    CXXFLAGS="-O3 -march=native" cmake -DBUILD_TESTING=ON -S . -B build --fresh
    cmake --build build --clean-first

And run tests:

    ctest --test-dir build --output-on-failure

## Examples

Please take a look at the examples:

* examples/isochrone
* examples/kepler
* examples/lotka-volterra

They will help you understand the intended logic of using the **collo**.

## Usage

The author currently suggests using **collo** as a git-submodule in your **CMake**-projects:

    ...
    add_subdirectory(external/collo)
    ...
    add_executable(executable_name ...)
    ...
    target_link_libraries(executable_name PRIVATE collo)
    ...

