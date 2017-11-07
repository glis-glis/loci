# -*- coding: utf-8 -*-
# loci
# Local cubic interpolations in up to 4 dimensions
#
# Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import ctypes as ct
import numpy as np
from . import CInterpolation

class Interpolation3D:
    """Interpolation in 3 dimensions. """

    def __init__(self, rx, ry, rz, f, f_x, f_y, f_z, f_xy, f_xz, f_yz, f_xyz):
        """Create interpolation in ranges `rx`, `ry`, `rz` for values 
        `fs`= f(`rx`, `ry`, `rz`) and derivative values `f_x` = df/dx, `f_y`,
        `f_z`, `f_xy`, `f_xz`, `f_yz`, and `f_xyz` = d^3f/(dx dy dz).
       
        `fs` is of length rx.len*ry.len. Its values are on every gridpoint.
        """
        self._lib                          = ct.CDLL("libloci.so")

        lib                          = self._lib
        lib.interpolation3D.restype  = CInterpolation
        lib.interpolate3D.restype    = ct.c_double
        lib.diff3D_x.restype         = ct.c_double
        lib.diff3D_y.restype         = ct.c_double
        lib.diff3D_z.restype         = ct.c_double
        lib.diff3D.restype           = ct.c_double
        lib.map3D.restype            = None
        lib.map3D_x.restype          = None
        lib.map3D_y.restype          = None
        lib.map3D_z.restype          = None
        lib.map3D_diff.restype       = None
        
        if callable(f): 
            #expects f AND f_x to be callable
            xsl     = np.r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx]
            ysl     = np.r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx]
            zsl     = np.r_[rz.x0:rz.x0 + rz.len*rz.dx:rz.dx]
            
            xs, ys, zs  = np.meshgrid(xsl, ysl, zsl, indexing='ij')
            xs          = xs.ravel()
            ys          = ys.ravel()
            zs          = zs.ravel()
            xyzs        = zip(xs, ys, zs)

            fs      = np.array([f(xi, yi, zi) for xi, yi, zi in xyzs])
            fs_x    = np.array([f_x(xi, yi, zi) for xi, yi, zi in xyzs])
            fs_y    = np.array([f_y(xi, yi, zi) for xi, yi, zi in xyzs])
            fs_z    = np.array([f_z(xi, yi, zi) for xi, yi, zi in xyzs])
            fs_xy   = np.array([f_xy(xi, yi, zi) for xi, yi, zi in xyzs])
            fs_xz   = np.array([f_xz(xi, yi, zi) for xi, yi, zi in xyzs])
            fs_yz   = np.array([f_yz(xi, yi, zi) for xi, yi, zi in xyzs])
            fs_xyz  = np.array([f_xyz(xi, yi, zi) for xi, yi, zi in xyzs])

        else:
            fs      = f
            fs_x    = f_x
            fs_y    = f_y
            fs_z    = f_z
            fs_xy   = f_xy
            fs_xz   = f_xz
            fs_yz   = f_yz
            fs_xyz  = f_xyz

        self._ip = lib.interpolation3D(ct.pointer(rx), ct.pointer(ry), 
                     ct.pointer(rz),
                     fs.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_x.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_y.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_z.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xy.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xz.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_yz.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xyz.ctypes.data_as(ct.POINTER(ct.c_double)))

    def interpolate(self, x, y, z):
        """Interpolate at point (`x`, `y`, `z`). """
        return self._lib.interpolate3D(ct.pointer(self._ip), 
                       ct.c_double(x), ct.c_double(y), ct.c_double(z))

    def diff_x(self, x, y, z):
        """Interpolate first order derivative in x at point (`x`, `y`, `z`). """
        return self._lib.diff3D_x(ct.pointer(self._ip), 
                       ct.c_double(x), ct.c_double(y), ct.c_double(z))

    def diff_y(self, x, y, z):
        """Interpolate first order derivative in y at point (`x`, `y`, `z`). """
        return self._lib.diff3D_y(ct.pointer(self._ip), 
                       ct.c_double(x), ct.c_double(y), ct.c_double(z))

    def diff_z(self, x, y, z):
        """Interpolate first order derivative in z at point (`x`, `y`, `z`). """
        return self._lib.diff3D_z(ct.pointer(self._ip), 
                       ct.c_double(x), ct.c_double(y), ct.c_double(z))

    def diff(self, edx, edy, edz, x, y, z):
        """Interpolate (`edx`, `edy`, `edz`)-order derivative in (x, y, z) at 
        point (`x`, `y`, `z`).
        """
        return self._lib.diff3D(ct.pointer(self._ip), ct.c_int(edx), 
                        ct.c_int(edy), ct.c_int(edz), 
                       ct.c_double(x), ct.c_double(y), ct.c_double(z))

    def map(self, xs, ys, zs):
        """Map interpolation onto the input-arrays `xs`, `ys`, `zs`
        and return result.
         
        The mapping is linear:
        result[i] = ip.interpolate(xs[i], ys[i], zs[i])
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map3D(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_x(self, xs, ys, zs):
        """Map first order derivative in x onto the input-arrays `xs`, `ys`, 
        `zs` and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map3D_x(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_y(self, xs, ys, zs):
        """Map first order derivative in y onto the input-arrays `xs`, `ys`, 
        `zs` and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map3D_y(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_z(self, xs, ys, zs):
        """Map first order derivative in z onto the input-arrays `xs`, `ys`, 
        `zs` and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map3D_z(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_diff(self, edx, edy, edz, xs, ys, zs):
        """Map (`edx`, `edy`, `edz`)-order derivative in (x, y, z) onto the
        input-arrays `xs`, `ys`, `zs` and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map3D_diff(ct.pointer(self._ip), ct.c_int(edx), ct.c_int(edy),
            ct.c_int(edz), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def __del__(self):
        self._lib.loci_free(ct.pointer(self._ip))
