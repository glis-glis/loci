# -*- coding: utf-8 -*-
# loci
# Local cubic interpolations in up to 4 dimensions
#
# Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from loci import Interpolation
from loci import Range
from numpy import *
from numpy.linalg import norm
from scipy import *

from helpers import tols, out

# Test functions and derivative.
funs    = [
    lambda x: 1.,
    lambda x: x,
    lambda x: 2*x + 1,
    lambda x: log(x*x + 1)
]

funs_x  = [
    lambda x: 0.,
    lambda x: 1.,
    lambda x: 2.,
    lambda x: 2*x/(x*x + 1)
]

def basic1D():
    """Test basic functionality:
    - create interpolation
    - interpolate a point
    - interpolate derivative
    - map interpolation
    - map derivation. 
    """

    out(0, "1D Basic Test")

    rg  = Range(1., 0.3, 20)                        # Interpolation range 

    xs0 = r_[rg.x0:rg.x0 + (rg.len-1)*rg.dx:rg.dx]  # Grid points
    xst = xs0 + rg.dx/3.                            # Test points

    xt  = rg.x0 + rg.dx*rg.len/2 + rg.dx/7          # Test value 

    for i in r_[0:len(funs)]:
        out(1, "Function %i"%(i+1))
        f   = funs[i]
        f_x = funs_x[i]

        # Exact values
        ft  = f(xt)
        ft_x= f_x(xt)

        # Create interpolatio
        ip  = Interpolation(rg, f, f_x)

        # Compare interpolated value to exact value
        assert abs(ip.interpolate(xt) - ft) < tols[1] + tols[1]*abs(ft),    \
               "interpolation(xt)"

        # Compare interpolated derivative to exact value
        assert abs(ip.diff_x(xt) - ft_x) < tols[1] + tols[1]*abs(ft_x),    \
               "diff(xt)"
        
        ###

        #  Map exact values on grid points 
        fs0     = array([f(xi) for xi in xs0])
        fs0_x   = array([f_x(xi) for xi in xs0])

        # Compare interpolated values to exact values on grid-points
        assert norm(fs0 - ip.map(xs0)) < tols[0] + tols[0]*norm(fs0),       \
               "Interpolation at xp=0"
        # Compare interpolated derivatives to exact values on grid-points
        assert norm(fs0_x - ip.map_x(xs0)) < tols[0] + tols[0]*norm(fs0_x), \
               "Diff_x at xp=0"

        ###

        # Map exact values on test points
        fst     = array([f(xi) for xi in xst])
        fst_x   = array([f_x(xi) for xi in xst])

        # Compare interpolated values to exact values on test-point
        assert norm(fst - ip.map(xst)) < tols[1] + tols[1]*norm(fst),       \
               "Interpolation in cell"

        # Compare interpolated derivatives to exact values on test-points
        assert norm(fst_x - ip.map_x(xst)) < tols[2] + tols[2]*norm(fst_x), \
               "Diff_x in cell"
    out(0, "Pass\n")

def exp1D():
    """Test higher order derivative on exponantial function. """

    out(0, "1D Exp Test");

    # See basic1D()
    rg  = Range(1., 0.3, 20)

    xs  = r_[rg.x0:rg.x0 + rg.len*rg.dx:rg.dx]
    xst = xs[0:-1] + rg.dx/3.

    xt  = rg.x0 + rg.dx*rg.len/2 + rg.dx/7
    ft  = exp(xt)

    fst = exp(xst)
    nfst= norm(fst)

    # Create interpolation
    ip  = Interpolation(rg, exp, exp)

    for i in r_[0:3]:
        out(1, "d^%if/dx^%i"%(i, i));

        ti  = tols[i+1]

        # Compare interpolated derivative to exact value
        assert abs(ip.diff(i, xt) - ft) < ti + ti*abs(ft), \
               "Diff in cell"

        # Compare interpolated derivatives to exact values on test-points
        assert norm(fst - ip.map_diff(i, xst)) < ti + ti*nfst, \
               "Diff in cell"

    out(0, "Pass\n")

def boundary1D():
    """Test out of bounds behaviour. """

    out(0, "1D Boundary Test")

    # See basic1D()
    big = 1e6
    fr  = 89./97

    rg  = Range(1., 0.3, 20)
    xs  = r_[rg.x0:rg.x0 + rg.len*rg.dx:rg.dx]

    # Points out of bounds 
    xts = [
        xs[-1]  + 1./big,
        xs[0]   - 1./big,
        xs[-1]  + rg.dx*fr,
        xs[0]   - rg.dx*fr,
        xs[-1]  + big,
        xs[0]   - big
    ]

    # Create interpolation
    ip  = Interpolation(rg, funs[-1], funs_x[-1])

    for xt in xts:
        # Test whether out of bounds == nan
        assert isnan(ip.interpolate(xt)), "Interpolation out of bounds"
        assert isnan(ip.diff_x(xt)), "Diff_x out of bounds"

    out(0, "Pass")
