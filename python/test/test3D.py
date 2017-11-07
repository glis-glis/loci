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

# Test functions and derivative
A   = 2
B   = 0.5
C   = 0.75

funs    = [
    lambda x, y, z: x + y + z,
    lambda x, y, z: x*y*z,
    lambda x, y, z: x**2*y**2*z**2,
    lambda x, y, z: log(A*x**2 + B*y**2 + + C*z**2 + 1)
]

funs_x  = [
    lambda x, y, z: 1.,
    lambda x, y, z: y*z,
    lambda x, y, z: 2*x*y**2*z**2,
    lambda x, y, z: 2*A*x/(A*x**2 + B*y**2 + C*z**2 + 1)
]

funs_y  = [
    lambda x, y, z: 1.,
    lambda x, y, z: x*z,
    lambda x, y, z: 2*x**2*y*z**2,
    lambda x, y, z: 2*B*y/(A*x**2 + B*y**2 + C*z**2 + 1)
]

funs_z  = [
    lambda x, y, z: 1.,
    lambda x, y, z: x*y,
    lambda x, y, z: 2*x**2*y**2*z,
    lambda x, y, z: 2*C*z/(A*x**2 + B*y**2 + C*z**2 + 1)
]

funs_xy = [
    lambda x, y, z: 0.,
    lambda x, y, z: z,
    lambda x, y, z: 4*x*y*z**2,
    lambda x, y, z: -4*A*B*x*y/(A*x**2 + B*y**2 + C*z**2 + 1)**2
]

funs_xz = [
    lambda x, y, z: 0.,
    lambda x, y, z: y,
    lambda x, y, z: 4*x*y**2*z,
    lambda x, y, z: -4*A*C*x*z/(A*x**2 + B*y**2 + C*z**2 + 1)**2
]

funs_yz = [
    lambda x, y, z: 0.,
    lambda x, y, z: x,
    lambda x, y, z: 4*x**2*y*z,
    lambda x, y, z: -4*B*C*y*z/(A*x**2 + B*y**2 + C*z**2 + 1)**2
]

funs_xyz = [
    lambda x, y, z: 0.,
    lambda x, y, z: 1,
    lambda x, y, z: 8*x*y*z,
    lambda x, y, z: 16*A*B*C*x*y*z/(A*x**2 + B*y**2 + C*z**2 + 1)**3
]

e3s =   [
    lambda x, y, z: exp(A*x + B*y + C*z),
    lambda x, y, z: A*exp(A*x + B*y + C*z),
    lambda x, y, z: B*exp(A*x + B*y + C*z),
    lambda x, y, z: C*exp(A*x + B*y + C*z),
    lambda x, y, z: A*B*exp(A*x + B*y + C*z),
    lambda x, y, z: A*C*exp(A*x + B*y + C*z),
    lambda x, y, z: B*C*exp(A*x + B*y + C*z),
    lambda x, y, z: A*B*C*exp(A*x + B*y + C*z)
]

def basic3D():
    """Test basic functionality:
    - create interpolation
    - interpolate a point
    - interpolate derivative
    - map interpolation
    - map derivation. 
    """

    out(0, "3D Basic Test")

    # Interpolation grid (rx*ry*rz) 
    rx  = Range(1., 0.3, 20)
    ry  = Range(-5., 0.2, 30)
    rz  = Range(2, 0.1, 10)

    # Grid points without last one
    xsl = r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx][0:-1]
    ysl = r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx][0:-1]
    zsl = r_[rz.x0:rz.x0 + rz.len*rz.dx:rz.dx][0:-1]

    xs0, ys0, zs0   = meshgrid(xsl, ysl, zsl, indexing='ij')
    xs0             = xs0.ravel()
    ys0             = ys0.ravel()
    zs0             = zs0.ravel()

    xyzs0           = zip(xs0, ys0, zs0)

    # Test points
    xst     = xs0 + rx.dx/3.
    yst     = ys0 + ry.dx/3.
    zst     = zs0 + rz.dx/3.
    xyzst   = zip(xst, yst, zst)

    # Test values 
    xt  = rx.x0 + rx.dx*rx.len/2 + rx.dx/7
    yt  = ry.x0 + ry.dx*ry.len/2 + ry.dx/7
    zt  = rz.x0 + rz.dx*rz.len/2 + rz.dx/7

    for i in r_[0:len(funs)]:
        out(1, "Function %i"%(i+1))
        f       = funs[i]
        f_x     = funs_x[i]
        f_y     = funs_y[i]
        f_z     = funs_z[i]
        f_xy    = funs_xy[i]
        f_xz    = funs_xz[i]
        f_yz    = funs_yz[i]
        f_xyz   = funs_xyz[i]

        # Exact values
        ft      = f(xt, yt, zt)
        ft_x    = f_x(xt, yt, zt)
        ft_y    = f_y(xt, yt, zt)
        ft_z    = f_z(xt, yt, zt)

        # Create interpolation
        ip  = Interpolation(rx, ry, rz, f, f_x, f_y, f_z, 
                            f_xy, f_xz, f_yz, f_xyz)

        # Compare interpolated value to exact value
        assert abs(ip.interpolate(xt, yt, zt) - ft) < tols[1]+tols[1]*abs(ft), \
               "interpolation(xt, yt, zt)"

        # Compare interpolated derivative to exact value
        assert abs(ip.diff_x(xt, yt, zt) - ft_x) < tols[1]+tols[1]*abs(ft_x), \
               "diff_x(xt, yt, zt)"

        assert abs(ip.diff_y(xt, yt, zt) - ft_y) < tols[1]+tols[1]*abs(ft_y), \
               "diff_y(xt, yt, zt)"

        assert abs(ip.diff_z(xt, yt, zt) - ft_z) < tols[1]+tols[1]*abs(ft_z), \
               "diff_z(xt, yt, zt)"

        ###

        # Map exact values on grid points (without last one,
        # as it is out of bounds) 
        fs0     = array([f(xi, yi, zi) for xi, yi, zi in xyzs0])
        fs0_x   = array([f_x(xi, yi, zi) for xi, yi, zi in xyzs0])
        fs0_y   = array([f_y(xi, yi, zi) for xi, yi, zi in xyzs0])
        fs0_z   = array([f_z(xi, yi, zi) for xi, yi, zi in xyzs0])

        fst     = array([f(xi, yi, zi) for xi, yi, zi in xyzst])
        fst_x   = array([f_x(xi, yi, zi) for xi, yi, zi in xyzst])
        fst_y   = array([f_y(xi, yi, zi) for xi, yi, zi in xyzst])
        fst_z   = array([f_z(xi, yi, zi) for xi, yi, zi in xyzst])

        # Compare interpolated values to exact values on grid-points
        assert norm(fs0 - ip.map(xs0, ys0, zs0)) < \
               tols[0] + tols[0]*norm(fs0), \
               "map(interpolation, xs0, ys0, zs0)" 

        # Compare interpolated derivatives to exact values on grid-points
        assert norm(fs0_x - ip.map_x(xs0, ys0, zs0)) < \
               tols[0] + tols[0]*norm(fs0_x), \
               "map(diff_x, xs0, ys0, zs0)"

        assert norm(fs0_y - ip.map_y(xs0, ys0, zs0)) < \
               tols[0] + tols[0]*norm(fs0_y), \
               "map(diff_y, xs0, ys0, zs0)"

        assert norm(fs0_z - ip.map_z(xs0, ys0, zs0)) < \
               tols[0] + tols[0]*norm(fs0_z), \
               "map(diff_z, xs0, ys0, zs0)"

        ###

        # Map exact values on test points 
        assert norm(fst - ip.map(xst, yst, zst)) < \
               tols[1] + tols[1]*norm(fst),\
               "map(interpolation, xst, yst, zst)"

        # Compare interpolated values to exact values on test-points
        assert norm(fst_x - ip.map_x(xst, yst, zst)) < \
               tols[1] + tols[1]*norm(fst_x),\
               "map(diff_x, xst, yst, zst)"

        # Compare interpolated derivatives to exact values on test-points
        assert norm(fst_y - ip.map_y(xst, yst, zst)) < \
               tols[1] + tols[1]*norm(fst_y),\
               "map(diff_y, xst, yst, zst)"

        assert norm(fst_z - ip.map_z(xst, yst, zst)) < \
               tols[1] + tols[1]*norm(fst_z),\
               "map(diff_z, xst, yst, zst)"

    out(0, "Pass\n")

def exp3D():
    """ Test higher order derivative on exponantial function. """

    out(0, "3D Exp Test");

    # See basic3D()
    rx  = Range(1., 0.3, 20)
    ry  = Range(-5., 0.2, 30)
    rz  = Range(2, 0.1, 10)

    xsl = r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx][0:-1]
    ysl = r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx][0:-1]
    zsl = r_[rz.x0:rz.x0 + rz.len*rz.dx:rz.dx][0:-1]

    xs0, ys0, zs0   = meshgrid(xsl, ysl, zsl, indexing='ij')
    xs0             = xs0.ravel()
    ys0             = ys0.ravel()
    zs0             = zs0.ravel()

    xyzs0           = zip(xs0, ys0, zs0)

    xst     = xs0 + rx.dx/3.
    yst     = ys0 + ry.dx/3.
    zst     = zs0 + rz.dx/3.
    xyzst   = zip(xst, yst, zst)

    xt  = rx.x0 + rx.dx*rx.len/2 + rx.dx/7
    yt  = ry.x0 + ry.dx*ry.len/2 + ry.dx/7
    zt  = rz.x0 + rz.dx*rz.len/2 + rz.dx/7

    ft  = e3s[0](xt, yt, zt)

    fs  = array([e3s[0](xi, yi, zi) for xi, yi, zi in xyzst])

    # Create interpolation
    ip  = Interpolation(rx, ry, rz, e3s[0], e3s[1], e3s[2], e3s[3], 
                        e3s[4], e3s[5], e3s[6], e3s[7])

    for i in r_[0:3]:
     for j in r_[0:3]:
      for k in r_[0:3]:
        out(1, "d^%if/dx^%idy^%idz^%i"%(i+j+k, i, j, k));

        ti  = tols[max(i, j, k)+1]

        Ai  = A**i
        Bj  = B**j
        Ck  = C**k

        fti = Ai*Bj*Ck*ft
        fsi = Ai*Bj*Ck*fs
        nfsi= norm(fsi)

        # Compare interpolated derivative to exact value
        assert abs(ip.diff(i, j, k, xt, yt, zt) - fti) < ti + ti*abs(fti), \
                "diff(xt, yt, zt)"

        # Compare interpolated derivatives to exact values on test-points
        assert norm(fsi - ip.map_diff(i, j, k, xst, yst, zst)) < ti + ti*nfsi, \
                "map(diff, i, j, k, xst, yst, zst)"

    out(0, "Pass\n")

def boundary3D():
    """ Test out of bounds behaviour. """

    out(0, "3D Boundary Test")

    big = 1e6
    fr  = 89./97

    # See basic3D()
    rx  = Range(1., 0.3, 20)
    ry  = Range(-5., 0.2, 30)
    rz  = Range(2, 0.1, 10)

    xs  = r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx]
    ys  = r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx]
    zs  = r_[rz.x0:rz.x0 + rz.len*rz.dx:rz.dx]

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

    zts = [
        zs[0]   + fr*rz.dx, # in bounds
        zs[-1]  + 1./big,
        zs[0]   - 1./big,
        zs[-1]  + fr*rz.dx,
        zs[0]   - fr*rz.dx,
        zs[-1]  + big,
        zs[0]   - big
    ]

    # Create interpolation
    ip  = Interpolation(rx, ry, rz, funs[-1], 
                        funs_x[-1], funs_y[-1], funs_z[-1],
                        funs_xy[-1], funs_xz[-1], funs_yz[-1], funs_xyz[-1])

    for i in r_[0:len(xts)]:
     for j in r_[0:len(yts)]:
      for k in r_[0:len(zts)]:
        if i+j+k > 0:
            xt  = xts[i]
            yt  = yts[j]
            zt  = zts[k]

            # Test whether out of bounds == nan
            assert isnan(ip.interpolate(xt, yt, zt)), \
                   "interpolate(xi, yj, zk) out of bounds"

            assert isnan(ip.diff_x(xt, yt, zt)), \
                   "diff_x(xi, yj, zk) out of bounds"

            assert isnan(ip.diff_y(xt, yt, zt)), \
                   "diff_y(xi, yj, zk) out of bounds"

            assert isnan(ip.diff_z(xt, yt, zt)), \
                   "diff_z(xi, yj, zk) out of bounds"

            assert isnan(ip.diff(2, 2, 2, xt, yt, zt)), \
                   "diff(2, 2, 2, xi, yj, zk) out of bounds"


    out(0, "Pass")
