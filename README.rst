======================
loci
======================
**Smooth Cubic Multivariate Local Interpolations**

|travis-ci|

loci is a shared library for interpolations in up to 4 dimensions. In order to
calculate the coefficients of the cubic polynom, only local values are used:
The data itself and all combinations of first-order derivatives, i.e. in 2D
f_x, f_y and f_xy. This is in contrast to splines, where the coefficients are
not calculated using derivatives, but non-local data, which can lead to
over-smoothing the result.

The scheme has been developed at the University of Geneva. It is based on
Lekien & Marsden 2005, with improvements by Daniel Pfenniger and implemented by
Andreas Füglistaler in C, Python and Julia.

.. |travis-ci| image:: https://api.travis-ci.com/glis-glis/loci.svg?branch=master
    :target: https://travis-ci.com/github/glis-glis/loci

Developement
============

The library has been developed using test-driven developement (TDD), i.e. a
test case is written **before** the actual function is implemented. That way,
one makes sure every test can fail, and that the actual function passes the
test.

The library has been tested with double precision. It is possible to compile
with single precision, by defining the variable `REAL=float`. This has however
not been tested yet.

The library can be linked both statically and dynamically. A python-wrapper is
available in `python/loci`, which is only working with double precision at this
moment.


Usage
============

Clone a copy of loci to your local computer and make the shared library:

.. code:: bash

    git clone https://github.com/AFueglistaler/loci.git
    cd loci
    make    

The shared library will be in `lib/libloci.so`, the header-files in `include/`,
the python module in `python/loci`. 

Testing
========

To test the statically linked test-files, type

.. code:: bash
    
    make test

To test the dynamically linked python-test-files, type

.. code:: bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/lib
    export PYTHONPATH=$PYTHONPATH:$(pwd)/python
    make pytest

Examples
========

The C-testfiles in `test/` and the python-testfiles in `python/test` are
commented and show the functionality of all functions.


Python Usage
============

In order to use loci in Python, you first need to set the `LD_LIBRARY_PATH` and
`PYTHONPATH` variables:

.. code:: bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path_to_loci/lib
    export PYTHONPATH=$PYTHONPATH:path_to_loci/python

The following code shows the basic usage of the library in two dimensions. The
usage in other dimensions is identical. 

.. code:: python

    #import numpy and scipy
    from numpy import *
    from scipy import *
    # import classes
    from loci import Interpolation, Range 

    # Define functions and derivative
    A=2.
    B=0.5
    def f(x, y):    return log(A*x**2 + B*y**2 + 1)
    def f_x(x, y):  return 2*A*x/(A*x**2 + B*y**2 + 1)
    def f_y(x, y):  return 2*B*y/(A*x**2 + B*y**2 + 1)
    def f_xy(x, y): return -4*A*B*x*y/(A*x**2 + B*y**2 + 1)**2
    
    # Define interpolation ranges
    rx  = Range(1., 0.1, 10)    #x0 =1., dx=0.1, lenght=10
    ry  = Range(-2., 0.5, 20)

    # Create interpolation
    ip  = Interpolation(rx, ry, f, f_x, f_y, f_xy)

    # Interpolate at a given point
    ip.interpolate(rx.x0 + 0.4, ry.x0 + 7.3)
    # Interpolate derivatives in x and y
    ip.diff_x(rx.x0 + 0.4, ry.x0 + 7.3)
    ip.diff_y(rx.x0 + 0.4, ry.x0 + 7.3)
    # Interpolate 2nd-order x and 3rd-order y derivative
    ip.diff(2, 3, rx.x0 + 0.4, ry.x0 + 7.3)

    # Interpolate out of bounds
    ip.interpolate(rx.x0 - 1, ry.x0 - 1)    #returns nan

    # create random points in ranges rx and ry
    N   = int(1e7)
    xs  = (rx.dx*rx.len)*rand(N) + rx.x0
    ys  = (ry.dx*ry.len)*rand(N) + ry.x0
    
    # Map interpolation on points
    ip.map(xs, ys)
    # Map derivativews in x and y on points
    ip.map_x(xs, ys)
    ip.map_y(xs, ys)
    # Map 2nd-order x and 3rd-order y derivative on points 
    ip.map_diff(2, 3, xs, ys)

Jupyter Notebooks
=================

There are two jupyter notebooks showing the usage of loci:

* `Introduction + Example <python/notebooks/Introduction%2BExample.ipynb>`_
* `Solid H2 Mass Fraction in the ISM <python/notebooks/Solid%20H2%20Mass%20Fraction%20in%20the%20ISM.ipynb>`_
