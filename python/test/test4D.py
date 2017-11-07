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

# Test functions
A   = 2
B   = 0.5
C   = 0.75
D   = 0.4

funs    = [
    lambda x, y, z, w: x + y + z + w,
    lambda x, y, z, w: x*y*z*w,
    lambda x, y, z, w: x**2*y**2*z**2*w**2,
    lambda x, y, z, w: log(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)
]

funs_x  = [
    lambda x, y, z, w: 1.,
    lambda x, y, z, w: y*z*w,
    lambda x, y, z, w: 2*x*y**2*z**2*w**2,
    lambda x, y, z, w: 2*A*x/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)
]

funs_y  = [
    lambda x, y, z, w: 1.,
    lambda x, y, z, w: x*z*w,
    lambda x, y, z, w: 2*x**2*y*z**2*w**2,
    lambda x, y, z, w: 2*B*y/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)
]

funs_z  = [
    lambda x, y, z, w: 1.,
    lambda x, y, z, w: x*y*w,
    lambda x, y, z, w: 2*x**2*y**2*z*w**2,
    lambda x, y, z, w: 2*C*z/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)
]

funs_w  = [
    lambda x, y, z, w: 1.,
    lambda x, y, z, w: x*y*z,
    lambda x, y, z, w: 2*x**2*y**2*z**2*w,
    lambda x, y, z, w: 2*D*w/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)
]

funs_xy = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: z*w,
    lambda x, y, z, w: 4*x*y*z**2*w**2,
    lambda x, y, z, w: -4*A*B*x*y/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**2
]

funs_xz = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: y*w,
    lambda x, y, z, w: 4*x*y**2*z*w**2,
    lambda x, y, z, w: -4*A*C*x*z/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**2
]

funs_xw = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: y*z,
    lambda x, y, z, w: 4*x*y**2*z**2*w,
    lambda x, y, z, w: -4*A*D*x*w/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**2
]

funs_yz = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: x*w,
    lambda x, y, z, w: 4*x**2*y*z*w**2,
    lambda x, y, z, w: -4*B*C*y*z/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**2
]

funs_yw = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: x*z,
    lambda x, y, z, w: 4*x**2*y*z**2*w,
    lambda x, y, z, w: -4*B*D*y*w/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**2
]

funs_zw = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: x*y,
    lambda x, y, z, w: 4*x**2*y**2*z*w,
    lambda x, y, z, w: -4*C*D*z*w/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**2
]

funs_xyz = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: w,
    lambda x, y, z, w: 8*x*y*z*w**2,
    lambda x, y, z, w: 16*A*B*C*x*y*z/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**3
]

funs_xyw = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: z,
    lambda x, y, z, w: 8*x*y*z**2*w,
    lambda x, y, z, w: 16*A*B*D*x*y*w/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**3
]

funs_xzw = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: y,
    lambda x, y, z, w: 8*x*y**2*z*w,
    lambda x, y, z, w: 16*A*C*D*x*z*w/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**3
]

funs_yzw = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: x,
    lambda x, y, z, w: 8*x**2*y*z*w,
    lambda x, y, z, w: 16*B*C*D*y*z*w/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**3
]

funs_xyzw = [
    lambda x, y, z, w: 0.,
    lambda x, y, z, w: 1.,
    lambda x, y, z, w: 16*x*y*z*w,
    lambda x, y, z, w: -96*A*B*C*D*x*y*z*w/(A*x**2 + B*y**2 + C*z**2 + D*w**2 + 1)**3
]

e4s =   [
    lambda x, y, z, w: exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: A*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: B*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: C*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: D*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: A*B*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: A*C*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: A*D*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: B*C*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: B*D*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: C*D*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: A*B*C*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: A*B*D*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: A*C*D*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: B*C*D*exp(A*x + B*y + C*z + D*w),
    lambda x, y, z, w: A*B*C*D*exp(A*x + B*y + C*z + D*w)
]

def basic4D():
    """ Test basic functionality:
    - create interpolation
    - interpolate a point
    - interpolate derivative
    - map interpolation
    - map derivation.
    """

    out(0, "4D Basic Test")

    # Interpolation grid (rx*ry*rz*rw)
    rx  = Range(1., 0.3, 20)
    ry  = Range(-5., 0.2, 30)
    rz  = Range(2, 0.1, 10)
    rw  = Range(-3, 0.4, 15)

    # Grid points without last one 
    xsl = r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx][0:-1]
    ysl = r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx][0:-1]
    zsl = r_[rz.x0:rz.x0 + rz.len*rz.dx:rz.dx][0:-1]
    wsl = r_[rw.x0:rw.x0 + rw.len*rw.dx:rw.dx][0:-1]

    xs0, ys0, zs0, ws0   = meshgrid(xsl, ysl, zsl, wsl, indexing='ij')
    xs0 = xs0.ravel()
    ys0 = ys0.ravel()
    zs0 = zs0.ravel()
    ws0 = ws0.ravel()

    xyzws0  = zip(xs0, ys0, zs0, ws0)

    # Test points
    xst     = xs0 + rx.dx/3.
    yst     = ys0 + ry.dx/3.
    zst     = zs0 + rz.dx/3.
    wst     = ws0 + rw.dx/3.
    xyzwst  = zip(xst, yst, zst, wst)

    # Test values
    xt  = rx.x0 + rx.dx*rx.len/2 + rx.dx/7
    yt  = ry.x0 + ry.dx*ry.len/2 + ry.dx/7
    zt  = rz.x0 + rz.dx*rz.len/2 + rz.dx/7
    wt  = rw.x0 + rw.dx*rw.len/2 + rw.dx/7

    for i in r_[0:len(funs)]:
        out(1, "Function %i"%(i+1))

        # Exact values
        f       = funs[i]
        f_x     = funs_x[i]
        f_y     = funs_y[i]
        f_z     = funs_z[i]
        f_w     = funs_w[i]
        f_xy    = funs_xy[i]
        f_xz    = funs_xz[i]
        f_xw    = funs_xw[i]
        f_yz    = funs_yz[i]
        f_yw    = funs_yw[i]
        f_zw    = funs_zw[i]
        f_xyz   = funs_xyz[i]
        f_xyw   = funs_xyw[i]
        f_xzw   = funs_xzw[i]
        f_yzw   = funs_yzw[i]
        f_xyzw  = funs_xyzw[i]

        ft      = f(xt, yt, zt, wt)
        ft_x    = f_x(xt, yt, zt, wt)
        ft_y    = f_y(xt, yt, zt, wt)
        ft_z    = f_z(xt, yt, zt, wt)
        ft_w    = f_w(xt, yt, zt, wt)
        
        # Create interpolation
        ip  = Interpolation(rx, ry, rz, rw, f, f_x, f_y, f_z, f_w,
                            f_xy, f_xz, f_xw, f_yz, f_yw, f_zw, 
                            f_xyz, f_xyw, f_xzw, f_yzw, f_xyzw)

        # Compare interpolated value to exact value 
        assert abs(ip.interpolate(xt, yt, zt, wt) - ft) < \
                    tols[1] + tols[1]*abs(ft), \
               "interpolation(xt, yt, zt, wt)"

        # Compare interpolated derivative to exact value
        assert abs(ip.diff_x(xt, yt, zt, wt) - ft_x) < \
                    tols[1] + tols[1]*abs(ft_x), \
               "interpolation(xt, yt, zt, wt)"

        assert abs(ip.diff_y(xt, yt, zt, wt) - ft_y) < \
                    tols[1] + tols[1]*abs(ft_y), \
               "interpolation(xt, yt, zt, wt)"

        assert abs(ip.diff_z(xt, yt, zt, wt) - ft_z) < \
                    tols[1] + tols[1]*abs(ft_z), \
               "interpolation(xt, yt, zt, wt)"

        assert abs(ip.diff_w(xt, yt, zt, wt) - ft_w) < \
                    tols[1] + tols[1]*abs(ft_w), \
               "interpolation(xt, yt, zt, wt)"

        ###
        
        # Map exact values on grid points (without last one,
        # as it is out of bounds)
        fs0     = array([f(xi, yi, zi, wi) for xi, yi, zi, wi in xyzws0])
        fs0_x   = array([f_x(xi, yi, zi, wi) for xi, yi, zi, wi in xyzws0])
        fs0_y   = array([f_y(xi, yi, zi, wi) for xi, yi, zi, wi in xyzws0])
        fs0_z   = array([f_z(xi, yi, zi, wi) for xi, yi, zi, wi in xyzws0])
        fs0_w   = array([f_w(xi, yi, zi, wi) for xi, yi, zi, wi in xyzws0])


        # Compare interpolated values to exact values on grid-points
        assert norm(fs0 - ip.map(xs0, ys0, zs0, ws0)) < \
               tols[0] + tols[0]*norm(fs0), \
               "map(interpolation, xs0, ys0, zs0, ws0)" 

        # Compare interpolated derivatives to exact values on grid-points
        assert norm(fs0_x - ip.map_x(xs0, ys0, zs0, ws0)) < \
               tols[0] + tols[0]*norm(fs0_x), \
               "map(diff_x, xs0, ys0, zs0, ws0)" 

        assert norm(fs0_y - ip.map_y(xs0, ys0, zs0, ws0)) < \
               tols[0] + tols[0]*norm(fs0_y), \
               "map(diff_y, xs0, ys0, zs0, ws0)" 

        assert norm(fs0_z - ip.map_z(xs0, ys0, zs0, ws0)) < \
               tols[0] + tols[0]*norm(fs0_z), \
               "map(diff_z, xs0, ys0, zs0, ws0)" 

        assert norm(fs0_w - ip.map_w(xs0, ys0, zs0, ws0)) < \
               tols[0] + tols[0]*norm(fs0_w), \
               "map(diff_w, xs0, ys0, zs0, ws0)" 

        ###

        # Map exact values on test points
        fst     = array([f(xi, yi, zi, wi) for xi, yi, zi, wi in xyzwst])
        fst_x   = array([f_x(xi, yi, zi, wi) for xi, yi, zi, wi in xyzwst])
        fst_y   = array([f_y(xi, yi, zi, wi) for xi, yi, zi, wi in xyzwst])
        fst_z   = array([f_z(xi, yi, zi, wi) for xi, yi, zi, wi in xyzwst])
        fst_w   = array([f_w(xi, yi, zi, wi) for xi, yi, zi, wi in xyzwst])

        # Compare interpolated values to exact values on test-points
        assert norm(fst - ip.map(xst, yst, zst, wst)) < \
               tols[1] + tols[1]*norm(fst), \
               "map(interpolation, xst, yst, zst, wst)" 

        # Compare interpolated derivatives to exact values on test-points
        assert norm(fst_x - ip.map_x(xst, yst, zst, wst)) < \
               tols[1] + tols[1]*norm(fst_x), \
               "map(diff_x, xst, yst, zst, wst)" 

        assert norm(fst_y - ip.map_y(xst, yst, zst, wst)) < \
               tols[1] + tols[1]*norm(fst_y), \
               "map(diff_y, xst, yst, zst, wst)" 

        assert norm(fst_z - ip.map_z(xst, yst, zst, wst)) < \
               tols[1] + tols[1]*norm(fst_z), \
               "map(diff_z, xst, yst, zst, wst)" 

        assert norm(fst_w - ip.map_w(xst, yst, zst, wst)) < \
               tols[1] + tols[1]*norm(fst_w), \
               "map(diff_w, xst, yst, zst, wst)" 

    out(0, "Pass\n")

def exp4D():
    """ Test higher order derivative on exponantial function. """
    out(0, "4D Exp Test");

    # See basic4D()
    rx  = Range(1., 0.3, 20)
    ry  = Range(-5., 0.2, 30)
    rz  = Range(2, 0.1, 10)
    rw  = Range(-3, 0.4, 15)

    xsl = r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx][0:-1]
    ysl = r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx][0:-1]
    zsl = r_[rz.x0:rz.x0 + rz.len*rz.dx:rz.dx][0:-1]
    wsl = r_[rw.x0:rw.x0 + rw.len*rw.dx:rw.dx][0:-1]

    xs0, ys0, zs0, ws0   = meshgrid(xsl, ysl, zsl, wsl, indexing='ij')
    xs0 = xs0.ravel()
    ys0 = ys0.ravel()
    zs0 = zs0.ravel()
    ws0 = ws0.ravel()

    xyzws0  = zip(xs0, ys0, zs0, ws0)

    xst     = xs0 + rx.dx/3.
    yst     = ys0 + ry.dx/3.
    zst     = zs0 + rz.dx/3.
    wst     = ws0 + rw.dx/3.
    xyzwst  = zip(xst, yst, zst, wst)

    xt  = rx.x0 + rx.dx*rx.len/2 + rx.dx/7
    yt  = ry.x0 + ry.dx*ry.len/2 + ry.dx/7
    zt  = rz.x0 + rz.dx*rz.len/2 + rz.dx/7
    wt  = rw.x0 + rw.dx*rw.len/2 + rw.dx/7

    ft  = e4s[0](xt, yt, zt, wt)

    fs  = array([e4s[0](xi, yi, zi, wi) for xi, yi, zi, wi in xyzwst])

    # Create interpolation
    ip  = Interpolation(rx, ry, rz, rw, e4s[0], e4s[1], e4s[2], e4s[3], 
                        e4s[4], e4s[5], e4s[6], e4s[7], e4s[8], e4s[9],
                        e4s[10], e4s[11], e4s[12], e4s[13], e4s[14], e4s[15])

    for i in r_[0:3]:
     for j in r_[0:3]:
      for k in r_[0:3]:
       for l in r_[0:3]:
        out(1, "d^%if/dx^%idy^%idz^%idw^%i"%(i+j+k, i, j, k, l));

        ti  = tols[max(i, j, k, l)+1]

        Ai  = A**i
        Bj  = B**j
        Ck  = C**k
        Dl  = D**l

        fti = Ai*Bj*Ck*Dl*ft
        fsi = Ai*Bj*Ck*Dl*fs
        nfsi= norm(fsi)

        # Compare interpolated derivative to exact value
        assert abs(ip.diff(i, j, k, l, xt, yt, zt, wt) - fti) < \
               ti + ti*abs(fti), \
               "diff(xt, yt, zt, wt)"

        # Compare interpolated derivatives to exact values on test-points
        assert norm(fsi - ip.map_diff(i, j, k, l, xst, yst, zst, wst)) < \
                ti + ti*nfsi, "map(diff, i, j, k, xst, yst, zst, wst)"

    out(0, "Pass\n")

def boundary4D():
    """ Test out of bounds behaviour. """

    out(0, "4D Boundary Test")

    big = 1e6
    fr  = 89./97

    # See basic4D()
    rx  = Range(1., 0.3, 20)
    ry  = Range(-5., 0.2, 30)
    rz  = Range(2, 0.1, 10)
    rw  = Range(-3, 0.4, 15)

    xs  = r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx]
    ys  = r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx]
    zs  = r_[rz.x0:rz.x0 + rz.len*rz.dx:rz.dx]
    ws  = r_[rw.x0:rw.x0 + rw.len*rw.dx:rw.dx]

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

    wts = [
        ws[0]   + fr*rw.dx, # in bounds
        ws[-1]  + 1./big,
        ws[0]   - 1./big,
        ws[-1]  + fr*rw.dx,
        ws[0]   - fr*rw.dx,
        ws[-1]  + big,
        ws[0]   - big
    ]

    # Create interpolation
    ip  = Interpolation(rx, ry, rz, rw, funs[-1], 
                        funs_x[-1], funs_y[-1], funs_z[-1], funs_w[-1],
                        funs_xy[-1], funs_xz[-1], funs_xw[-1], funs_yz[-1], 
                        funs_yw[-1], funs_zw[-1], funs_xyz[-1], funs_xyw[-1], 
                        funs_xzw[-1], funs_yzw[-1], funs_xyzw[-1])

    for i in r_[0:len(xts)]:
     for j in r_[0:len(yts)]:
      for k in r_[0:len(zts)]:
       for l in r_[0:len(wts)]:
        if i+j+k+l > 0:
            xt  = xts[i]
            yt  = yts[j]
            zt  = zts[k]
            wt  = wts[l]

            # Test whether out of bounds == nan
            assert isnan(ip.interpolate(xt, yt, zt, wt)), \
                   "interpolate(xi, yj, zk, wl) out of bounds"

            assert isnan(ip.diff_x(xt, yt, zt, wt)), \
                   "diff_x(xi, yj, zk, wl) out of bounds"

            assert isnan(ip.diff_y(xt, yt, zt, wt)), \
                   "diff_y(xi, yj, zk, wl) out of bounds"

            assert isnan(ip.diff_z(xt, yt, zt, wt)), \
                   "diff_z(xi, yj, zk, wl) out of bounds"

            assert isnan(ip.diff_w(xt, yt, zt, wt)), \
                   "diff_w(xi, yj, zk, wl) out of bounds"

            assert isnan(ip.diff(2, 2, 2, 2, xt, yt, zt, wt)), \
                   "diff(xi, yj, zk, wl) out of bounds"

    out(0, "Pass")
