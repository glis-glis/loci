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
A   = 2.
B   = 0.5

funs    = [
    lambda x, y: x + y,
    lambda x, y: x*y,
    lambda x, y: x**2*y**2,
    lambda x, y: log(A*x**2 + B*y**2 + 1)
]

funs_x  = [
    lambda x, y: 1.,
    lambda x, y: y,
    lambda x, y: 2*x*y**2,
    lambda x, y: 2*A*x/(A*x**2 + B*y**2 + 1)
]

funs_y  = [
    lambda x, y: 1.,
    lambda x, y: x,
    lambda x, y: 2*x**2*y,
    lambda x, y: 2*B*y/(A*x**2 + B*y**2 + 1)
]

funs_xy = [
    lambda x, y: 0.,
    lambda x, y: 1.,
    lambda x, y: 4*x*y,
    lambda x, y: -4*A*B*x*y/(A*x**2 + B*y**2 + 1)**2
]

e2s =   [
    lambda x, y: exp(A*x+B*y),
    lambda x, y: A*exp(A*x+B*y),
    lambda x, y: B*exp(A*x+B*y),
    lambda x, y: A*B*exp(A*x+B*y)
]

def basic2D():
    """Test basic functionality:
    - create interpolation
    - interpolate a point
    - interpolate derivative
    - map interpolation
    - map derivation. 
    """

    out(0, "2D Basic Test")

    # Interpolation grid (rx * ry) 
    rx  = Range(1., 0.3, 20)
    ry  = Range(-5., 0.2, 30)

    # Grid points without last one
    xsl = r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx][0:-1]
    ysl = r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx][0:-1]

    xs0, ys0    = meshgrid(xsl, ysl, indexing='ij')
    xs0         = xs0.ravel()
    ys0         = ys0.ravel()
    xys0        = zip(xs0, ys0)

    # Test points
    xst = xs0 + rx.dx/3.
    yst = ys0 + ry.dx/3.
    xyst= zip(xst, yst)

    # Test values
    xt  = rx.x0 + rx.dx*rx.len/2 + rx.dx/7
    yt  = ry.x0 + ry.dx*ry.len/2 + ry.dx/7

    for i in r_[0:len(funs)]:
        out(1, "Function %i"%(i+1))

        f       = funs[i]
        f_x     = funs_x[i]
        f_y     = funs_y[i]
        f_xy    = funs_xy[i]

        # Exact values
        ft      = f(xt, yt)
        ft_x    = f_x(xt, yt)
        ft_y    = f_y(xt, yt)
        ft_xy   = f_xy(xt, yt)

        # Create interpolation
        ip  = Interpolation(rx, ry, f, f_x, f_y, f_xy)

        # Compare interpolated value to exact value
        assert abs(ip.interpolate(xt, yt) - ft) < tols[1] + tols[1]*abs(ft), \
               "interpolation(xt, yt)"

        # Compare interpolated derivative to exact value
        assert abs(ip.diff_x(xt, yt) - ft_x) < tols[1] + tols[1]*abs(ft_x), \
               "diff_x(xt, yt)"

        assert abs(ip.diff_y(xt, yt) - ft_y) < tols[1] + tols[1]*abs(ft_y), \
               "diff_y(xt, yt)"

        ###

        #  Map exact values on grid points (without last one,
        # as it is out of bounds) 
        fs0     = array([f(xi, yi) for xi, yi in xys0])
        fs0_x   = array([f_x(xi, yi) for xi, yi in xys0])
        fs0_y   = array([f_y(xi, yi) for xi, yi in xys0])

        # Compare interpolated values to exact values on grid-points
        assert norm(fs0 - ip.map(xs0, ys0)) < tols[0] + tols[0]*norm(fs0), \
               "map(interpolation, xs0, ys0)" 

        # Compare interpolated derivatives to exact values on grid-points
        assert norm(fs0_x-ip.map_x(xs0, ys0)) < tols[0]+tols[0]*norm(fs0_x), \
               "map(diff_x, xs0, ys0)"

        assert norm(fs0_y-ip.map_y(xs0, ys0)) < tols[0]+tols[0]*norm(fs0_x), \
               "map(diff_y, xs0, ys0)"

        ###

        # Map exact values on test points
        fst     = array([f(xi, yi) for xi, yi in xyst])
        fst_x   = array([f_x(xi, yi) for xi, yi in xyst])
        fst_y   = array([f_y(xi, yi) for xi, yi in xyst])

        # Compare interpolated values to exact values on test-points
        assert norm(fst - ip.map(xst, yst)) < tols[1] + tols[1]*norm(fst),\
               "map(interpolation, xst, yst)"

        # Compare interpolated derivatives to exact values on test-points
        assert norm(fst_x-ip.map_x(xst, yst)) < tols[2]+tols[2]*norm(fst_x), \
               "map(diff_x, xst, yst)"

        assert norm(fst_y-ip.map_y(xst, yst)) < tols[2]+tols[2]*norm(fst_y), \
               "map(diff_y, xst, yst)"

    out(0, "Pass\n")

def exp2D():
    """Test higher order derivative on exponantial function. """

    out(0, "2D Exp Test");

    # See basic2D()
    rx  = Range(1., 0.3, 20)
    ry  = Range(-5., 0.2, 30)

    xsl = r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx][0:-1]
    ysl = r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx][0:-1]

    xs0, ys0    = meshgrid(xsl, ysl, indexing='ij')
    xs0         = xs0.ravel()
    ys0         = ys0.ravel()
    xys0        = zip(xs0, ys0)

    xst = xs0 + rx.dx/3.
    yst = ys0 + ry.dx/3.
    xyst= zip(xst, yst)

    xt  = rx.x0 + rx.dx*rx.len/2 + rx.dx/7
    yt  = ry.x0 + ry.dx*ry.len/2 + ry.dx/7

    ft  = e2s[0](xt, yt)

    fs  = array([e2s[0](xi, yi) for xi, yi in xyst])

    # Create interpolation
    ip  = Interpolation(rx, ry, e2s[0], e2s[1], e2s[2], e2s[3])

    for i in r_[0:3]:
        for j in r_[0:3]:
            out(1, "d^%if/dx^%idy^%i"%(i+j, i, j));

            ti  = tols[max(i, j)+1]

            Ai  = A**i
            Bj  = B**j

            fti = Ai*Bj*ft
            fsi = Ai*Bj*fs
            nfsi= norm(fsi)

            # Compare interpolated derivative to exact value
            assert abs(ip.diff(i, j, xt, yt) - fti) < ti + ti*abs(fti), \
                    "diff(xt, yt)"

            # Compare interpolated derivatives to exact values on 
            # test-points 
            assert norm(fsi - ip.map_diff(i, j, xst, yst)) < ti + ti*nfsi, \
                    "map(diff, i, j, xst, yst)"

    out(0, "Pass\n")

def boundary2D():
    """Test out of bounds behaviour. """

    out(0, "2D Boundary Test")

    big = 1e6
    fr  = 89./97

    # See basic2D() 
    rx  = Range(1., 0.3, 20)
    ry  = Range(-5., 0.2, 30)

    xs  = r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx]
    ys  = r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx]

    # Points out of bounds 
    xts = [
        xs[0]   + fr*rx.dx, # in bounds
        xs[-1]  + 1./big,
        xs[0]   - 1./big,
        xs[-1]  + fr*rx.dx,
        xs[0]   - fr*rx.dx,
        xs[-1]  + big,
        xs[0]   - big
    ]

    yts = [
        ys[0]   + fr*ry.dx, # in bounds
        ys[-1]  + 1./big,
        ys[0]   - 1./big,
        ys[-1]  + fr*ry.dx,
        ys[0]   - fr*ry.dx,
        ys[-1]  + big,
        ys[0]   - big
    ]

    # Create interpolation
    ip  = Interpolation(rx, ry, funs[-1], funs_x[-1], funs_y[-1], funs_xy[-1])

    for i in r_[0:len(xts)]:
        for j in r_[0:len(yts)]:
            if i+j > 0:
                xt  = xts[i]
                yt  = yts[j]

                # Test whether out of bounds == nan
                assert isnan(ip.interpolate(xt, yt)), \
                       "interpolate(xi, yj) out of bounds"
                assert isnan(ip.diff_x(xt, yt)), "diff_x(xi, yj) out of bounds"
                assert isnan(ip.diff_y(xt, yt)), "diff_y(xi, yj) out of bounds"
                assert isnan(ip.diff(2, 2, xt, yt)), \
                       "diff(2, 2, xi, yj) out of bounds"

    out(0, "Pass")
