# collo

This is the header-only library for collocation integrator for ODEs (ordinary differential equations).

The integrator offers an adjustable order, single-step, implicit method for solving ODEs : *y' = f(t, y)*.
Here *y* is a n-dimensional vector that is represented by **Eigen::Matrix**.

You can find details about collocation integrators in these books:
* Preston C. Hammer and J W Hollingsworth. “Trapezoidal methods of approximating solutions of differential equations”. B: Mathematics of Computation 9 (1955), c. 92—96.
* Ernst Hairer, Christian Lubich and Gerhard Wanner. “Geometric numerical integration illustrated by the Störmer–Verlet method”. B: Acta Numerica 12 (2003), c. 399—450.

## Getting started

You can define *f(t, y)* as a function object with two arguments: time *t* and vector *y*.
Use **auto** as a type of second argument or use the helper alias **collo::state_vector_t<some_float_type, *n*>**.

Include **<collo/lobatto.hpp>** and then you can create an integrator with the function **make_Lobatto**. Example:

    `auto lobatto = make_Lobatto<some_float_type, n, s, type_of_f>(y_0, t_0, h);`

or

    `auto lobatto = make_Lobatto<some_float_type, n, s>(y_0, t_0, h, f);`

Parameter *s* determines how many intermediate points are used for the calculation of result on every step. This parameter impacts accuracy of integrator $O(h^{2s-2})$.
Arguments *y_0* and *t_0* are initial vector and time. Argument *h* is a time step.

### Methods of the integrator class.

Use `do_step()` to integrate a single step with a fixed step size. <br>

* `steps()` returns the number of calculated steps.
* `time()` returns the current time *t = t_0 + h* · `steps()`.
* `iternum()` returns the number of iterations of the main implicit method loop on the last step.
* `state()` returns the current state *y*.
* `force()` can be used to get access to *f* and change it if needed.
* `poly(some_float_type t)` returns the intermediate interpolated vector *y* with the argument `t` being between 0 and 1. Other values of argument will give extrapolation with collocation polynomial.
* `poly_node(size_t i)` returns the vector *y* at collocation node with number i.

## Compilation

This project is dependent on the library `Eigen` (version 3.4+). You can get it from [official page](https://eigen.tuxfamily.org/index.php?title=Main_Page).

    `git clone https://github.com/shvak/collo.git`
    `cd collo`
    `mkdir build`
    `CXXFLAGS="-O3 -march=native" cmake -DBUILD_TESTING=OFF -S . -B build`
    `cmake --build build`

Note: You can improve the performance of the integrator if you compile your code with the option '-march=native'. This feature is provided by the 'Eigen' library and modern CPU vectorization.

Additionally, if you want to use built-in tests, you will need the library `doctest`. Then:

    `CXXFLAGS="-O3 -march=native" cmake -DBUILD_TESTING=ON -S . -B build --fresh`
    `cmake --build build --clean-first`

And run tests:

    `ctest --test-dir build --output-on-failure`

## Examples

Please take a look at the examples:

    `build/examples/isochrone/isochrone`
    `build/examples/kepler/kepler`
    `build/examples/lotka-volterra/lotka-volterra`

They will help you understand the intended logic of using the **collo**.

## Usage

The author suggests using the **collo** as a git-submodule at the moment.
