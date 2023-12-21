# collo

This is the repo for collocation integrator for ODEs (ordinary differential equations).

The integrator offers an adjustable order, single-step, implicit method for solving ODEs : *y' = f(t, y)*.
Here *y* is a n-dimensional vector that is represented by **Eigen::Matrix**.

You can find details about collocation integrators in these books:
* Preston C. Hammer and J W Hollingsworth. “Trapezoidal methods of approximating solutions of differential equations”. B: Mathematics of Computation 9 (1955), c. 92—96.
* Ernst Hairer, Christian Lubich and Gerhard Wanner. “Geometric numerical integration illustrated by the Störmer–Verlet method”. B: Acta Numerica 12 (2003), c. 399—450.

## Getting started

You can define *f(t, y)* as a function object with two arguments: time *t* and vector *y*.
For better performance, use **auto** as a type of second argument. Or you can use the helper alias **collo::state_vector_t<some_float_type, *n*>**.

Include **<collo/lobatto.hpp>** or **<collo/gauss.hpp>** and then you can create an integrator with the functions **make_Lobatto** or **make_Gauss**. Example:

`auto lobatto = make_Lobatto<some_float_type, n, s, type_of_f>(y_0, t_0, h);`

Parameter *s* determines how many intermediate points are used for the calculation of result on every step. This parameter impacts accuracy of integrator $O(h^{2s-2})$.
Arguments *y_0* and *t_0* are initial vector and time. Argument *h* is a time step.

### Methods of the integrator class. 
Use `do_step()` to integrate a single step with a fixed step size. <br>

* `steps()` returns the number of calculated steps.
* `time()` returns the current timepoint.
* `iternum()` returns the number of iterations of the main implicit method loop on the last step.
* `state()` returns the current state *y*.
* `force()` can be used to get access to *f* and change it if needed.
* `poly(some_float_type t)` returns the intermediate interpolated vector *y* with the argument `t` being between 0 and 1. Other values of argument will give extrapolation with collocation polynomial.
*  `poly_node(size_t i)` returns the vector *y* at collocation node with number i.




## Compilation

Note: You can improve the performance of the integrator if you compile your code with the option '-march=native'. This feature is provided by the 'Eigen' library and modern CPU vectorization. Like this:

`git clone https://github.com/shvak/collo.git`<br>
 `cd collo`<br>
 `mkdir build`<br>
 `CXXFLAGS="-O3 -march=native" cmake -S . -B build`<br>
 `cmake --build build`<br>

Run simple tests:

` ctest --test-dir build --output-on-failure`<br>

## Examples

`build/examples/isochrone/isochrone`<br>
`build/examples/kepler/kepler`<br>
`build/examples/lotka-volterra/lotka-volterra`<br>
