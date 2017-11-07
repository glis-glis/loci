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

class Interpolation4D:
    """Interpolation in 4 dimensions. """

    def __init__(self, rx, ry, rz, rw, f, f_x, f_y, f_z, f_w, 
                 f_xy, f_xz, f_xw, f_yz, f_yw, f_zw, 
                 f_xyz, f_xyw, f_xzw, f_yzw, f_xyzw):
        """Create interpolation in ranges `rx`, `ry`, `rz`, `rw` for values 
        `fs`= f(`rx`, `ry`, `rz`) and derivative values `f_x` = df/dx, `f_y`,
        `f_z`, `f_w`, `f_xy`, `f_xz`, `f_xw`, `f_yz`, `f_yw`, `f_zw`, `f_xyz`,
        `f_xyz`, `f_xyw`,`f_xzw`, `f_yzw`, `f_xyzw`= d^4f/(dx dy dz dw). 

        `fs` is of length rx.len*ry.len*rz.len*rw.len. 
        Its values are on every gridpoint.
        """
        self._lib                          = ct.CDLL("libloci.so")

        lib                          = self._lib
        lib.interpolation4D.restype  = CInterpolation
        lib.interpolate4D.restype    = ct.c_double
        lib.diff4D_x.restype         = ct.c_double
        lib.diff4D_y.restype         = ct.c_double
        lib.diff4D_z.restype         = ct.c_double
        lib.diff4D_w.restype         = ct.c_double
        lib.diff4D.restype           = ct.c_double
        lib.map4D.restype            = None
        lib.map4D_x.restype          = None
        lib.map4D_y.restype          = None
        lib.map4D_z.restype          = None
        lib.map4D_w.restype          = None
        lib.map4D_diff.restype       = None
        
        if callable(f): 
            #expects f AND f_x to be callable
            xsl     = np.r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx]
            ysl     = np.r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx]
            zsl     = np.r_[rz.x0:rz.x0 + rz.len*rz.dx:rz.dx]
            wsl     = np.r_[rw.x0:rw.x0 + rw.len*rw.dx:rw.dx]
            
            xs, ys, zs, ws  = np.meshgrid(xsl, ysl, zsl, wsl, indexing='ij')
            xs      = xs.ravel()
            ys      = ys.ravel()
            zs      = zs.ravel()
            ws      = ws.ravel()
            xyzws   = zip(xs, ys, zs, ws)

            fs      = np.array([f(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_x    = np.array([f_x(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_y    = np.array([f_y(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_z    = np.array([f_z(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_w    = np.array([f_w(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_xy   = np.array([f_xy(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_xz   = np.array([f_xz(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_xw   = np.array([f_xw(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_yz   = np.array([f_yz(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_yw   = np.array([f_yw(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_zw   = np.array([f_zw(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_xyz  = np.array([f_xyz(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_xyw  = np.array([f_xyw(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_xzw  = np.array([f_xzw(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_yzw  = np.array([f_yzw(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])
            fs_xyzw = np.array([f_xyzw(xi, yi, zi, wi) 
                                for xi, yi, zi, wi in xyzws])

        else:
            fs      = f
            fs_x    = f_x
            fs_y    = f_y
            fs_z    = f_z
            fs_w    = f_w
            fs_xy   = f_xy
            fs_xz   = f_xz
            fs_xw   = f_xw
            fs_yz   = f_yz
            fs_yw   = f_yw
            fs_zw   = f_zw
            fs_xyz  = f_xyz
            fs_xyw  = f_xyw
            fs_xzw  = f_xzw
            fs_yzw  = f_yzw
            fs_xyzw = f_xyzw

        self._ip = lib.interpolation4D(ct.pointer(rx), ct.pointer(ry), 
                     ct.pointer(rz), ct.pointer(rw),
                     fs.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_x.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_y.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_z.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_w.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xy.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xz.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xw.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_yz.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_yw.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_zw.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xyz.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xyw.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xzw.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_yzw.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xyzw.ctypes.data_as(ct.POINTER(ct.c_double)))

    def interpolate(self, x, y, z, w):
        """Interpolate at point (`x`, `y`, `z`, `w`). """
        return self._lib.interpolate4D(ct.pointer(self._ip), 
                 ct.c_double(x), ct.c_double(y), ct.c_double(z), ct.c_double(w))

    def diff_x(self, x, y, z, w):
        """Interpolate first order derivative in x
        at point (`x`, `y`, `z`, `w`). 
        """
        return self._lib.diff4D_x(ct.pointer(self._ip), 
                 ct.c_double(x), ct.c_double(y), ct.c_double(z), ct.c_double(w))

    def diff_y(self, x, y, z, w):
        """Interpolate first order derivative in y
        at point (`x`, `y`, `z`, `w`). 
        """
        return self._lib.diff4D_y(ct.pointer(self._ip), 
                 ct.c_double(x), ct.c_double(y), ct.c_double(z), ct.c_double(w))

    def diff_z(self, x, y, z, w):
        """Interpolate first order derivative in z
        at point (`x`, `y`, `z`, `w`). 
        """
        return self._lib.diff4D_z(ct.pointer(self._ip), 
                 ct.c_double(x), ct.c_double(y), ct.c_double(z), ct.c_double(w))

    def diff_w(self, x, y, z, w):
        """Interpolate first order derivative in w
        at point (`x`, `y`, `z`, `w`). 
        """
        return self._lib.diff4D_w(ct.pointer(self._ip), 
                 ct.c_double(x), ct.c_double(y), ct.c_double(z), ct.c_double(w))

    def diff(self, edx, edy, edz, edw, x, y, z, w):
        """Interpolate (`edx`, `edy`, `edz`, `edw`)-order derivative in 
        (x, y, z, w) at point (`x`, `y`, `z`, `w`).
        """
        return self._lib.diff4D(ct.pointer(self._ip), ct.c_int(edx), 
                ct.c_int(edy), ct.c_int(edz), ct.c_int(edw), 
                ct.c_double(x), ct.c_double(y), ct.c_double(z), ct.c_double(w))

    def map(self, xs, ys, zs, ws):
        """Map interpolation onto the input-arrays `xs`, `ys`, `zs`, `ws`
        and return result.
         
        The mapping is linear:
        result[i] = ip.interpolate(xs[i], ys[i], zs[i], ws[i])
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map4D(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ws.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_x(self, xs, ys, zs, ws):
        """Map first order derivative in x onto the input-arrays `xs`, `ys`, 
        `zs`, `ws` and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map4D_x(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ws.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_y(self, xs, ys, zs, ws):
        """Map first order derivative in y onto the input-arrays `xs`, `ys`, 
        `zs`, `ws` and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map4D_y(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ws.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_z(self, xs, ys, zs, ws):
        """Map first order derivative in z onto the input-arrays `xs`, `ys`, 
        `zs`, `ws` and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map4D_z(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ws.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_w(self, xs, ys, zs, ws):
        """Map first order derivative in w onto the input-arrays `xs`, `ys`, 
        `zs`, `ws` and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map4D_w(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ws.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_diff(self, edx, edy, edz, edw, xs, ys, zs, ws):
        """Map (`edx`, `edy`, `edz`, `edw`)-order derivative in (x, y, z) onto 
        the input-arrays `xs`, `ys`, `zs`, `ws` and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map4D_diff(ct.pointer(self._ip), ct.c_int(edx), ct.c_int(edy),
            ct.c_int(edz), ct.c_int(edw), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            zs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ws.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def __del__(self):
        self._lib.loci_free(ct.pointer(self._ip))
