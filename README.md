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
Method **state()** can be used to get current *y*.<br>
Method **force()** can be used to get access to *f* and change it if needed.
Method **poly()** can be used to get intermediate vector *y* with argument between 0 and 1. Other values of argument will give extrapolation with collocation polynom.
Method **poly_node()** can be used to get vector *y* at collocation nodes.

Please see examples for better understanding.

Note: You can rich full power of integrator if compile your code and provided examples with option '-march=native'. This feature provided by 'Eigen' library and modern cpu's vectorization.
