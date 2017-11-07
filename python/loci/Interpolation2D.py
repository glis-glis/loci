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

class Interpolation2D:
    """Interpolation in 2 dimensions. """ 

    def __init__(self, rx, ry, f, f_x, f_y, f_xy):
        """Create interpolation in ranges `rx`, `ry` for values 
         `fs`= f(`rx`, `ry`) and derivative values `f_x` = df/dx, `f_y` and 
         `f_xy` = d^2f/(dx dy). 
         
         `fs` is of length rx.len*ry.len. Its values are on every gridpoint.
        """
        self._lib                          = ct.CDLL("libloci.so")

        lib                          = self._lib
        lib.interpolation2D.restype  = CInterpolation
        lib.interpolate2D.restype    = ct.c_double
        lib.diff2D_x.restype         = ct.c_double
        lib.diff2D_y.restype         = ct.c_double
        lib.diff2D.restype           = ct.c_double
        lib.map2D.restype            = None
        lib.map2D_x.restype          = None
        lib.map2D_y.restype          = None
        lib.map2D_diff.restype       = None
        
        if callable(f): 
            #expects all f* to be callable
            xsl     = np.r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx]
            ysl     = np.r_[ry.x0:ry.x0 + ry.len*ry.dx:ry.dx]
            
            xs, ys  = np.meshgrid(xsl, ysl, indexing='ij')
            xs      = xs.ravel()
            ys      = ys.ravel()
            xys     = zip(xs, ys)

            fs      = np.array([f(xi, yi) for xi, yi in xys])
            fs_x    = np.array([f_x(xi, yi) for xi, yi in xys])
            fs_y    = np.array([f_y(xi, yi) for xi, yi in xys])
            fs_xy   = np.array([f_xy(xi, yi) for xi, yi in xys])

        else:
            fs      = f
            fs_x    = f_x
            fs_y    = f_y
            fs_xy   = f_xy

        self._ip = lib.interpolation2D(ct.pointer(rx), ct.pointer(ry),
                     fs.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_x.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_y.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_xy.ctypes.data_as(ct.POINTER(ct.c_double)))

    def interpolate(self, x, y):
        """Interpolate at point (`x`, `y`). """
        return self._lib.interpolate2D(ct.pointer(self._ip), 
                       ct.c_double(x), ct.c_double(y))

    def diff_x(self, x, y):
        """Interpolate first order derivative in x at point (`x`, `y`). """
        return self._lib.diff2D_x(ct.pointer(self._ip), 
                       ct.c_double(x), ct.c_double(y))

    def diff_y(self, x, y):
        """Interpolate first order derivative in y at point (`x`, `y`). """
        return self._lib.diff2D_y(ct.pointer(self._ip), 
                       ct.c_double(x), ct.c_double(y))

    def diff(self, edx, edy, x, y):
        """Interpolate (`edx`, `edy`)-order derivative in (x, y) at 
        point (`x`, `y`).
        """
        return self._lib.diff2D(ct.pointer(self._ip), ct.c_int(edx), 
                        ct.c_int(edy), ct.c_double(x), ct.c_double(y)) 

    def map(self, xs, ys):
        """Map interpolation onto the input-arrays `xs`, `ys` 
        and return result.
        
        The mapping is linear:
        result[i] = ip.interpolate(xs[i], ys[i])
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map2D(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_x(self, xs, ys):
        """Map first order derivative in x onto the input-arrays `xs`, `ys` 
        and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map2D_x(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_y(self, xs, ys):
        """Map first order derivative in y onto the input-arrays `xs`, `ys` 
        and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map2D_y(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out


    def map_diff(self, edx, edy, xs, ys):
        """Map (`edx`, `edy`)-order derivative in (x, y) onto the 
        input-arrays `xs`, `ys` and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map2D_diff(ct.pointer(self._ip), ct.c_int(edx), ct.c_int(edy),
            ct.c_int(l), xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            ys.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def __del__(self):
        self._lib.loci_free(ct.pointer(self._ip))
