# collo

Collocation integrator for ODEs: *y' = f(t, y)*.

Here *y* is a vector with size *n* and represented by **Eigen::Matrix**.
You can define *f(t, y)* as function object with two arguments: a time *t* and a vector *y*.
For better performance use **auto** as type of second argument. Or you can use helper alias **collo::state_vector_t<some_float_type, *n*>**.

Include **<collo/lobatto.hpp>** or **<collo/gauss.hpp>** and than you can create integrator with functions **make_Lobatto** or **make_Gauss**. Example:

> **auto lobatto = make_Lobatto<some_float_type, *n*, *s*, type_of_*f*>(*y_0*, *t_0*, *h*);**

Parameter *s* define how many intermediate points used for calculation of result on every step. Its value impact accuracy of integrator.
Arguments *y_0* and *t_0* are initial vector and time. Argument *h* is a time step.

Methods of integrator do a job. Use **do_step()** for calculations and use **steps()**, **time()**, **iternum()** and **state()** to get current state of integrator.<br>
Method **steps()** give number of calculated steps.<br>
Method **time()** allow to get current timepoint.<br>
Method **iternum()** give number of iterations of main cycle of last calculated step.<br>
Method **state()** must be used to get current *y*.

Please see examples for better understanding.

Note: C++20 standard code, that can be compiled by gcc only now, gcc version 11.0+. Clang wants too many 'typename's and doesn't want some consteval functions.

Note: You can rich full power of integrator if compile your code and provided examples with option '-march=native'. This feature provided by 'Eigen' library and modern cpu's vectorization.
